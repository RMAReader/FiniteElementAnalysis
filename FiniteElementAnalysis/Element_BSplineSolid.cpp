#include "stdafx.h"
#include "Element_BSplineSolid.h"
#include "BSplineSolid.h"
#include "bspline_utils.h"
#include "iostream"
#include "gauss_legendre.h"
#include <omp.h>
#include "slu_ddefs.h"
#include <vector>

Element_BSplineSolid::Element_BSplineSolid(BSplineSolid &newsolid, double* newD, double newdensity)
{
	solid = newsolid;
	density = newdensity;
	for (int i = 0; i < 9; i++){
		D[i] = newD[i];
	}

}


Element_BSplineSolid::~Element_BSplineSolid()
{
	//try{
	//	Destroy_CompCol_Matrix(&M);
	//	Destroy_CompCol_Matrix(&K);
	//}
	//catch (std::exception e){
	//	std::cout << "failed to destroy M and K" << std::endl;
	//}
}



SuperMatrix  &Element_BSplineSolid::getM(){

	Element_BSplineSolid::M_build(solid, density, &M);

	return M;
}


double Element_BSplineSolid::get_mass(){
	return 1.0;
}


SuperMatrix &Element_BSplineSolid::getK(){

	//transform elasticities into form used for calc
	double D1[27] =
	{ D[0], D[6], D[7], D[6], D[1], 0, D[7], 0, D[2],
	D[1], D[6], 0, D[6], D[3], D[8], 0, D[8], D[4],
	D[2], 0, D[7], 0, D[4], D[8], D[7], D[8], D[5] };

	Element_BSplineSolid::K_build(solid, D1, &K);

	return K;
}



int Element_BSplineSolid::buildBasisFunctions(double knot[], int lknot, int p, double raw_abiscas[], double  raw_weights[], int pts, double db[], double db_1[],double weights[]){
	
	double* values = new double[p];
	double a, b, abisca, weight;

	//iterate over each knot span
	int index = 0;
	for (int span = 0; span < lknot - 2 * p-1; span++){

		a = 0.5 * (knot[span + p + 1] + knot[span + p]);
		b = 0.5 * (knot[span + p + 1] - knot[span + p]);

		//iterate over each integration point within knot span
		for (int ip = 0; ip < pts; ip++){
			
			abisca = a + b * raw_abiscas[ip];
			weight = b * raw_weights[ip];

			//iterate over each point whose basis funtion is non-zero in knot span
			for (int cp = span; cp <= span + p; cp++){

				bspline::deboor_values(cp, p, knot, lknot, abisca, values);
				db[index] = values[0];
				db_1[index] = values[1];
				weights[index] = weight;
				index++;

			}

		}

	}

	return 0;
}


int Element_BSplineSolid::getGaussPoints(double knot[], int lknot, int p, int pts, double raw_abiscas[],double  raw_weights[],double abiscas[], double weights[],int lresult){
	
	if (lresult != (lknot - 2 * p)*pts) return -1;

	for (int i = 0; i < lknot - 2*p; i++){
		
		double a = 0.5 * (knot[i+p + 1] + knot[i+p]);
		double b = 0.5 * (knot[i+p + 1] - knot[i+p]);
		
		for (int j = 0; j < pts; j++){
			abiscas[i * pts + j] = a + b * raw_abiscas[j];
			weights[i * pts + j] = b * raw_weights[j];
			//std::cout << i * pts + j << ": " <<abiscas[i * pts + j] << ", " << weights[i * pts + j] << "\n";
		}
	}
	return 0;
}



int Element_BSplineSolid::M_build(BSplineSolid &solid, double density, SuperMatrix* M){

	Element_BSplineSolid::initialise_SparseCRSMatrix(solid, M);

	int ptx = 3 * solid.p();
	int pty = 3 * solid.q();
	int ptz = 3 * solid.r();

	double *raw_abiscas_x = new double[ptx];
	double *raw_abiscas_y = new double[pty];
	double *raw_abiscas_z = new double[ptz];
	double *raw_weights_x = new double[ptx];
	double *raw_weights_y = new double[pty];
	double *raw_weights_z = new double[ptz];
	gauss_legendre_tbl(ptx, raw_abiscas_x, raw_weights_x, 1E-15);
	gauss_legendre_tbl(pty, raw_abiscas_y, raw_weights_y, 1E-15);
	gauss_legendre_tbl(ptz, raw_abiscas_z, raw_weights_z, 1E-15);

	int g_points_x = solid.spansx() * solid.pspanx() * ptx;			//number of spans * control point which affect each span * integration points per span
	int g_points_y = solid.spansy() * solid.pspany() * pty;			//number of spans * control point which affect each span * integration points per span
	int g_points_z = solid.spansz() * solid.pspanz() * ptz;			//number of spans * control point which affect each span * integration points per span

	double* dbx = new double[g_points_x];		
	double* db_1x = new double[g_points_x];	
	double* weights_x = new double[g_points_x];	
	double* dby = new double[g_points_y];		
	double* db_1y = new double[g_points_y];	
	double* weights_y = new double[g_points_y];	
	double* dbz = new double[g_points_z];		
	double* db_1z = new double[g_points_z];	
	double* weights_z = new double[g_points_z];	

	Element_BSplineSolid::buildBasisFunctions(solid.knotx, solid.lknotx(), solid.p(), raw_abiscas_x, raw_weights_x, ptx, dbx, db_1x, weights_x);
	Element_BSplineSolid::buildBasisFunctions(solid.knoty, solid.lknoty(), solid.q(), raw_abiscas_y, raw_weights_y, pty, dby, db_1y, weights_y);
	Element_BSplineSolid::buildBasisFunctions(solid.knotz, solid.lknotz(), solid.r(), raw_abiscas_z, raw_weights_z, ptz, dbz, db_1z, weights_z);

	omp_set_num_threads(8);
	#pragma omp parallel for
	for (int a = 0; a < solid.spansx(); a++){
		for (int b = 0; b < solid.spansy(); b++){
			for (int c = 0; c < solid.spansz(); c++){

				//temporary array for dense results of single knot span
				double *values = new double[solid.pspantot() * solid.pspantot()];

				Element_BSplineSolid::M_buildSpan(solid,
					density,
					a, ptx, raw_weights_x, dbx, db_1x,
					b, pty, raw_weights_y, dby, db_1y,
					c, ptz, raw_weights_z, dbz, db_1z,
					values);

				#pragma omp critical
				{
					Element_BSplineSolid::M_addSpan(solid, a, b, c, values, M);
				}
					
				delete[] values;

			}
		}
		std::cout << "a: " << a << "\n";
	}

	delete[] raw_abiscas_x;
	delete[] raw_abiscas_y;
	delete[] raw_abiscas_z;
	delete[] raw_weights_x;
	delete[] raw_weights_y;
	delete[] raw_weights_z;
	delete[] dbx;
	delete[] db_1x;
	delete[] weights_x;
	delete[] dby;
	delete[] db_1y;
	delete[] weights_y;
	delete[] dbz;
	delete[] db_1z;
	delete[] weights_z;
	
	return 0;
}





int Element_BSplineSolid::M_buildSpan(BSplineSolid &solid, double density,
	int span_x, int ptx, double raw_weights_x[], double db_x[], double db_x1[],
	int span_y, int pty, double raw_weights_y[], double db_y[], double db_y1[],
	int span_z, int ptz, double raw_weights_z[], double db_z[], double db_z1[],
	double values[]){

	double* JacobianDet = new double[ptx*pty*ptz];
	double val0, val1, val2, val3, val4, val5, val6, val7; // , val8;
	int Jcounter,Pcounter;

	//offsets for db_ indices
	int span_xoffset = span_x * ptx * solid.pspanx();
	int span_yoffset = span_y * pty * solid.pspany();
	int span_zoffset = span_z * ptz * solid.pspanz();
	int dbx_indOffset;
	int dby_indOffset;
	int dbz_indOffset;

	double *weights_x = new double[ptx];
	double *weights_y = new double[pty];
	double *weights_z = new double[ptz];
	ScaleWeights(solid.knotx[span_x + solid.p()], solid.knotx[span_x + solid.p() + 1], ptx, raw_weights_x, weights_x);
	ScaleWeights(solid.knoty[span_y + solid.q()], solid.knoty[span_y + solid.q() + 1], pty, raw_weights_y, weights_y);
	ScaleWeights(solid.knotz[span_z + solid.r()], solid.knotz[span_z + solid.r() + 1], ptz, raw_weights_z, weights_z);

	//find Jacobian Determinants
	Element_BSplineSolid::buildJacobianDet(solid, 
		span_x, ptx, db_x, db_x1, 
		span_y, pty, db_y, db_y1, 
		span_z, ptz, db_z, db_z1, JacobianDet);

	//for (int dummy = 0; dummy < ptx*pty*ptz; dummy++){
	//	std::cout << "JacobianDet[" << dummy << "] = " << JacobianDet[dummy] <<"\n";
	//}


	//arrays for Kahan sum, used for increased flaoting point precision
	for (int i = 0; i < solid.pspantot() * solid.pspantot(); i++){
		values[i] = 0;
	}

	dbx_indOffset = span_xoffset;
	Jcounter = 0;
	for (int ipx = 0; ipx < ptx; ipx++, dbx_indOffset += solid.pspanx()){
				
		val0 = density * weights_x[ipx];
		
		dby_indOffset = span_yoffset;

		for (int ipy = 0; ipy < pty; ipy++, dby_indOffset += solid.pspany()){
			
			val1 = val0 * weights_y[ipy];
			
			dbz_indOffset = span_zoffset;

			for (int ipz = 0; ipz < ptz; ipz++, dbz_indOffset += solid.pspanz()){
				
				
				val2 = val1 * weights_z[ipz] * JacobianDet[Jcounter];
				Jcounter++;

				//cycle over rows of M Matrix
				Pcounter = 0;
				for (int i = dbx_indOffset; i <= dbx_indOffset+solid.p(); i++){

					val3 = val2 * db_x[i];

					for (int j = dby_indOffset; j <= dby_indOffset + solid.q(); j++){

						val4 = val3 * db_y[j];

						for (int k = dbz_indOffset; k <= dbz_indOffset + solid.r(); k++){

							val5 = val4 * db_z[k];

							//cycle over cols of M matrix, only covering upper triangular part
							for (int l = dbx_indOffset; l <= dbx_indOffset + solid.p(); l++){
								
								if (false) continue;

								val6 = val5 * db_x[l];

								for (int m = dby_indOffset; m <= dby_indOffset + solid.q(); m++){
									
									if (false) continue;

									val7 = val6 * db_y[m];
									
									for (int n = dbz_indOffset; n <= dbz_indOffset + solid.r(); n++){
										
										if (false) continue;

										values[Pcounter] += val7 * db_z[n];

										Pcounter++;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	
	delete[] JacobianDet;
	delete[] weights_x;
	delete[] weights_y;
	delete[] weights_z;

	return 0;
	
}


int Element_BSplineSolid::M_addSpan(BSplineSolid &solid, int span_x, int span_y, int span_z, double values[], SuperMatrix* M){
	//add results to M
	int Mindex1, Mindex2, Mindex3, row;
	int index = 0;

	NCformat *Mstore = (NCformat*) M->Store;
	double* nzval = (double*) Mstore->nzval;

	//for rows, cycle over points in knot span (span_x, span_y, span_z)
	for (int i = 0; i <= solid.p(); i++){
		for (int j = 0; j <= solid.q(); j++){
			for (int k = 0; k <= solid.r(); k++){

				row = 3 * ((i + span_x) * solid.coffsetx() + (j + span_y) * solid.coffsety() + (k + span_z * solid.coffsetz()));

				//for columns, find position of columns in compressed row
				int stepz = 3;
				int stepy = (min(solid.nz(), span_z + k + solid.r() + 1) - max(0, span_z + k - solid.r())) * stepz;
				int stepx = (min(solid.ny(), span_y + j + solid.q() + 1) - max(0, span_y + j - solid.q())) * stepy;

				for (int l = min(span_x, solid.p() - i); l <= min(span_x, solid.p() - i) + solid.p(); l++){
					for (int m = min(span_y, solid.q() - j); m <= min(span_y, solid.q() - j) + solid.q(); m++){
						for (int n = min(span_z, solid.r() - k); n <= min(span_z, solid.r() - k) + solid.r(); n++){

							Mindex1 = Mstore->colptr[row] + l * stepx + m * stepy + n * stepz;
							Mindex2 = Mstore->colptr[row + 1] + l * stepx + m * stepy + n * stepz + 1;
							Mindex3 = Mstore->colptr[row + 2] + l * stepx + m * stepy + n * stepz + 2;

							nzval[Mindex1] += values[index];
							nzval[Mindex2] += values[index];
							nzval[Mindex3] += values[index];

							index++;
						}
					}
				}
			}
		}
	}

	return 0;
}



/*

function to create array of JacobianDeterminants for all integration points in a given knot span

*/
int Element_BSplineSolid::buildJacobianDet(BSplineSolid &solid, 
	int span_x, int ptx, double db_x[], double db_x1[],
	int span_y, int pty, double db_y[], double db_y1[],
	int span_z, int ptz, double db_z[], double db_z1[],
	double JacobianDet[]){

	double Jacobian[9], J1[3];
	
	int p0 ,p1, p2, p3;
	int dx0, dy0, dz0, dx1, dy1, dz1;
	int counter;

	p0 = solid.coffsetx() * span_x;
	counter = 0;

	dx0 = span_x * ptx * solid.pspanx();
	for (int ipx = 0; ipx < ptx; ipx++, dx0 += solid.pspanx()){
		dy0 = span_y * pty * solid.pspany();
		for (int ipy = 0; ipy < pty; ipy++, dy0 += solid.pspany()){
			dz0 = span_z * ptz * solid.pspanz();
			for (int ipz = 0; ipz < ptz; ipz++, dz0 += solid.pspanz()){

				for (int dummy = 0; dummy < 9; dummy++){
					Jacobian[dummy] = 0;
				}
				p1 = p0;
				dx1 = dx0;

				for (int i = 0; i <= solid.p(); i++, dx1 += 1, p1 += solid.coffsetx()){

					p2 = p1 + solid.coffsety() * span_y;
					dy1 = dy0;

					for (int j = 0; j <= solid.q(); j++, dy1 += 1, p2 += solid.coffsety()){

						J1[0] = db_x1[dx1] * db_y[dy1];
						J1[1] = db_x[dx1] * db_y1[dy1];
						J1[2] = db_x[dx1] * db_y[dy1];

						p3 = p2 + solid.coffsetz() *span_z;
						dz1 = dz0;

						for (int k = 0; k <= solid.r(); k++, dz1 += 1, p3 += solid.coffsetz()){

							Jacobian[0] += J1[0] * db_z[dz1] * solid.cx[p3];
							Jacobian[1] += J1[0] * db_z[dz1] * solid.cy[p3];
							Jacobian[2] += J1[0] * db_z[dz1] * solid.cz[p3];
							Jacobian[3] += J1[1] * db_z[dz1] * solid.cx[p3];
							Jacobian[4] += J1[1] * db_z[dz1] * solid.cy[p3];
							Jacobian[5] += J1[1] * db_z[dz1] * solid.cz[p3];
							Jacobian[6] += J1[2] * db_z1[dz1] * solid.cx[p3];
							Jacobian[7] += J1[2] * db_z1[dz1] * solid.cy[p3];
							Jacobian[8] += J1[2] * db_z1[dz1] * solid.cz[p3];
							
						}
					}
				}

				JacobianDet[counter] = Element_BSplineSolid::Determinant(Jacobian);
				counter++;

				
			}
		}
	}

	return 0;
}

void Element_BSplineSolid::ScaleWeights(double lower, double upper, int n, double raw_weights[], double weights[]){
	double b = 0.5 * (upper - lower);
	for (int i = 0; i < n; i++){
		weights[i] = b * raw_weights[i];
	}
}

double Element_BSplineSolid::Determinant(double M[]){
	return M[0] * (M[4] * M[8] - M[5] * M[7]) - M[1] * (M[3] * M[8] - M[5] * M[6]) + M[2] * (M[3] * M[7] - M[4] * M[6]);
}

void Element_BSplineSolid::DetermnantInverse(double M[], double inv[], double &det){
	double a = (M[4] * M[8] - M[5] * M[7]);
	double b = (M[5] * M[6] - M[3] * M[8]);
	double c = (M[3] * M[7] - M[4] * M[6]);
	det = M[0] * a + M[1] * b + M[2] * c;

	inv[0] = a / det;
	inv[1] = (M[2] * M[7] - M[1] * M[8]) / det;
	inv[2] = (M[1] * M[5] - M[2] * M[4]) / det;
	inv[3] = b / det;
	inv[4] = (M[0] * M[8] - M[2] * M[6]) / det;
	inv[5] = (M[2] * M[3] - M[0] * M[5]) / det;
	inv[6] = c / det;
	inv[7] = (M[1] * M[6] - M[0] * M[7]) / det;
	inv[8] = (M[0] * M[4] - M[1] * M[3]) / det;
}

void Element_BSplineSolid::Product(double M[], double v[], double Mv[]){
	for (int i = 0, k = 0; i < 3; i++){
		Mv[i] = 0;
		for (int j = 0; j < 3; j++, k++){
			Mv[i] += M[k] * v[j];
		}
	}
}


void Element_BSplineSolid::initialise_SparseCRSMatrix(BSplineSolid &solid, SuperMatrix* matrix){

	Dtype_t dtype = SLU_D;
	Mtype_t mtype = SLU_GE;
	Stype_t stype = SLU_NC;
	int nrow = solid.npoints() * 3;
	int ncol = solid.npoints() * 3;

	int nnz = 9 * (solid.nx()*(2 * solid.p() + 1) - solid.p()*(solid.p() + 1))
					* (solid.ny()*(2 * solid.q() + 1) - solid.q()*(solid.q() + 1))
					* (solid.nz()*(2 * solid.r() + 1) - solid.r()*(solid.r() + 1));

	double* nzval = new double[nnz];
	int* rowind = new int[nnz];
	int* colptr = new int[nrow + 1];
	colptr[nrow] = nnz;

	int row = 0;
	int valCounter = 0;
	for (int i = 0; i < solid.nx(); i++){
		for (int j = 0; j < solid.ny(); j++){
			for (int k = 0; k < solid.nz(); k++){
				for (int u = 0; u < 3; u++, row++){

					colptr[row] = valCounter;
					
					for (int l = max(0, i - solid.p()); l <= min(solid.nx()-1, i + solid.p()); l++){
						for (int m = max(0, j - solid.q()); m <= min(solid.ny()-1, j + solid.q()); m++){
							for (int n = max(0, k - solid.r()); n <= min(solid.nz()-1, k + solid.r()); n++){
								for (int v = 0; v < 3; v++){
									
									nzval[valCounter] = 0;
									rowind[valCounter] = (l*solid.coffsetx() + m*solid.coffsety() + n*solid.coffsetz()) * 3 + v;
									valCounter++;
									
								}
							}
						}
					}
				}
			}
		}
	}

	dCreate_CompCol_Matrix(matrix, nrow, ncol, nnz, nzval, rowind, colptr, stype, dtype, mtype);
}



int Element_BSplineSolid::K_buildSpan(BSplineSolid &solid, double D[],
	int span_x, int ptx, double raw_weights_x[], double db_x[], double db_x1[],
	int span_y, int pty, double raw_weights_y[], double db_y[], double db_y1[],
	int span_z, int ptz, double raw_weights_z[], double db_z[], double db_z1[],
	double values[]){

	double dNu[3],dNu0[3];
	double *dNx = new double[4 * solid.pspanx()*solid.pspany()*solid.pspanz()];

	int Jcounter, Pcounter;

	//offsets for db_ indices
	int span_xoffset = span_x * ptx * solid.pspanx();
	int span_yoffset = span_y * pty * solid.pspany();
	int span_zoffset = span_z * ptz * solid.pspanz();
	int dbx_indOffset;
	int dby_indOffset;
	int dbz_indOffset;

	double *weights_x = new double[ptx];
	double *weights_y = new double[pty];
	double *weights_z = new double[ptz];
	ScaleWeights(solid.knotx[span_x + solid.p()], solid.knotx[span_x + solid.p() + 1], ptx, raw_weights_x, weights_x);
	ScaleWeights(solid.knoty[span_y + solid.q()], solid.knoty[span_y + solid.q() + 1], pty, raw_weights_y, weights_y);
	ScaleWeights(solid.knotz[span_z + solid.r()], solid.knotz[span_z + solid.r() + 1], ptz, raw_weights_z, weights_z);

	//arrays for Kahan sum, used for increased flaoting point precision
	for (int i = 0; i < 9 * solid.pspantot() * solid.pspantot(); i++){
		values[i] = 0;
	}

	//indx for finding dNx given position in D - 3 represents 0
	int index[27] = { 0, 1, 2, 1, 0, 3, 2, 3, 0,
		1, 0, 3, 0, 1, 2, 3, 2, 1,
		2, 3, 0, 3, 2, 1, 0, 1, 2 };

	double BD[27];   //array is ordered so can be multiplied in blocks of 3 by {dx,dy,dz} to give values
	// = {dx.D11, dy.D44, dz.D55,   dy.D44, dx.D12,      0,   dz.D55,      0, dx.D13,
	//	   dy.D21, dx.D44,      0,   dx.D44, dy.D22, dz.D66,        0, dz.D66, dy.D23,
	//     dz.D31,      0, dx.D55,        0, dz.D32, dy.D66,   dx.D55, dy.D66, dz.D33 }


	dbx_indOffset = span_xoffset;
	int p0 = solid.coffsetx() * span_x;
	int p1, p2, p3;
	Jcounter = 0;
	double h, h0;

	//cycle over integration points
	for (int ipx = 0; ipx < ptx; ipx++, dbx_indOffset += solid.pspanx()){

		dby_indOffset = span_yoffset;

		for (int ipy = 0; ipy < pty; ipy++, dby_indOffset += solid.pspany()){

			h0 = weights_x[ipx] * weights_y[ipy];
			dbz_indOffset = span_zoffset;

			for (int ipz = 0; ipz < ptz; ipz++, dbz_indOffset += solid.pspanz()){

				Jcounter++;
				//cycle over rows of M Matrix, building array of dNx(i,j,k) for (i,j,k) = 0 to (p,q,r)
				Pcounter = 0;


				//build Jacobian
				double Jacobian[9];
				double J1[3];
				double JacobianInverse[9];
				double JacobianDet;

				int dx1 = dbx_indOffset;
				int dy1, dz1;
				p1 = p0;

				for (int dummy = 0; dummy < 9; dummy++){
					Jacobian[dummy] = 0;
				}


				for (int i = dbx_indOffset; i <= dbx_indOffset + solid.p(); i++, p1 +=solid.coffsetx()){
					
					p2 = p1 + solid.coffsety() * span_y;

					for (int j = dby_indOffset; j <= dby_indOffset + solid.q(); j++, p2 +=solid.coffsety()){

						J1[0] = db_x1[i] * db_y[j];
						J1[1] = db_x[i] * db_y1[j];
						J1[2] = db_x[i] * db_y[j];

						p3 = p2 + solid.coffsetz() *span_z;

						for (int k = dbz_indOffset; k <= dbz_indOffset + solid.r(); k++, p3 +=solid.coffsetz()){

							Jacobian[0] += J1[0] * db_z[k] * solid.cx[p3];
							Jacobian[1] += J1[0] * db_z[k] * solid.cy[p3];
							Jacobian[2] += J1[0] * db_z[k] * solid.cz[p3];
							Jacobian[3] += J1[1] * db_z[k] * solid.cx[p3];
							Jacobian[4] += J1[1] * db_z[k] * solid.cy[p3];
							Jacobian[5] += J1[1] * db_z[k] * solid.cz[p3];
							Jacobian[6] += J1[2] * db_z1[k] * solid.cx[p3];
							Jacobian[7] += J1[2] * db_z1[k] * solid.cy[p3];
							Jacobian[8] += J1[2] * db_z1[k] * solid.cz[p3];

						}
					}
				}

				DetermnantInverse(Jacobian, JacobianInverse, JacobianDet);

				Jcounter = 0;
				for (int i = dbx_indOffset; i <= dbx_indOffset + solid.p(); i++){

					for (int j = dby_indOffset; j <= dby_indOffset + solid.q(); j++){

						dNu0[0] = db_x1[i] * db_y[j];
						dNu0[1] = db_x[i] * db_y1[j];
						dNu0[2] = db_x[i] * db_y[j];

						for (int k = dbz_indOffset; k <= dbz_indOffset + solid.r(); k++){

							dNu[0] = dNu0[0] * db_z[k];
							dNu[1] = dNu0[1] * db_z[k];
							dNu[2] = dNu0[2] * db_z1[k];

							Product(JacobianInverse, dNu, &dNx[Jcounter]);

							//double dNx1 = dNx[Jcounter];
							//double dNx2 = dNx[Jcounter+1];
							//double dNx3 = dNx[Jcounter+2];
							dNx[Jcounter + 3] = 0;
							Jcounter += 4;
						}
					}
				}
							
				int row = 0;
				int col = 0;
				Pcounter = 0;

				h = h0 * weights_z[ipz] * JacobianDet;
				for (int i = 0; i <= solid.p(); i++){

					for (int j = 0; j <= solid.q(); j++){

						for (int k = 0; k <= solid.r(); k++, row+=4){

							for (int u = 0; u < 27; u++){
								BD[u] = h * D[u] * dNx[row + index[u]];
							}
							
							//for (int dummy = 0; dummy < 27; dummy++){
							//	std::cout << "BD["<<dummy<<"] = "<<BD[dummy]<<"\n";
							//}


							//cycle over cols of M matrix, only covering upper triangular part
							col = 0;
							for (int l = 0; l <= solid.p(); l++){

								if (false) continue;
								
								for (int m = 0; m <= solid.q(); m++){

									if (false) continue;
									
									for (int n = 0; n <= solid.r(); n++, col+=4){

										if (false) continue;
										double tmp4 = 0;

										for (int u = Pcounter, w = 0; u < Pcounter + 9; u++){
											for (int v = 0; v < 3; v++, w++){

												//double tmp1 = dNx[col + v];
												//double tmp2 = BD[w];
												//double tmp3 = tmp1 * tmp2;
												//tmp4 += tmp3;

												values[u] += BD[w] * dNx[col + v];
											}
										}

										//Product(BD, &B[col], &values[Pcounter]);

										Pcounter += 9;
									}
								}
							}
						}
					}
				}
			}
		}
	}

	delete[] weights_x;
	delete[] weights_y;
	delete[] weights_z;
	delete[] dNx;

	return 0;

}




int Element_BSplineSolid::K_build(BSplineSolid &solid, double D[], SuperMatrix *K){

	Element_BSplineSolid::initialise_SparseCRSMatrix(solid, K);

	int ptx = 3 * solid.p();
	int pty = 3 * solid.q();
	int ptz = 3 * solid.r();

	double *raw_abiscas_x = new double[ptx];
	double *raw_abiscas_y = new double[pty];
	double *raw_abiscas_z = new double[ptz];
	double *raw_weights_x = new double[ptx];
	double *raw_weights_y = new double[pty];
	double *raw_weights_z = new double[ptz];
	gauss_legendre_tbl(ptx, raw_abiscas_x, raw_weights_x, 1E-15);
	gauss_legendre_tbl(pty, raw_abiscas_y, raw_weights_y, 1E-15);
	gauss_legendre_tbl(ptz, raw_abiscas_z, raw_weights_z, 1E-15);

	int g_points_x = solid.spansx() * solid.pspanx() * ptx;			//number of spans * control point which affect each span * integration points per span
	int g_points_y = solid.spansy() * solid.pspany() * pty;			//number of spans * control point which affect each span * integration points per span
	int g_points_z = solid.spansz() * solid.pspanz() * ptz;			//number of spans * control point which affect each span * integration points per span

	double* dbx = new double[g_points_x];
	double* db_1x = new double[g_points_x];
	double* weights_x = new double[g_points_x];
	double* dby = new double[g_points_y];
	double* db_1y = new double[g_points_y];
	double* weights_y = new double[g_points_y];
	double* dbz = new double[g_points_z];
	double* db_1z = new double[g_points_z];
	double* weights_z = new double[g_points_z];

	Element_BSplineSolid::buildBasisFunctions(solid.knotx, solid.lknotx(), solid.p(), raw_abiscas_x, raw_weights_x, ptx, dbx, db_1x, weights_x);
	Element_BSplineSolid::buildBasisFunctions(solid.knoty, solid.lknoty(), solid.q(), raw_abiscas_y, raw_weights_y, pty, dby, db_1y, weights_y);
	Element_BSplineSolid::buildBasisFunctions(solid.knotz, solid.lknotz(), solid.r(), raw_abiscas_z, raw_weights_z, ptz, dbz, db_1z, weights_z);

	omp_set_num_threads(8);
#pragma omp parallel for
	for (int a = 0; a < solid.spansx(); a++){
		for (int b = 0; b < solid.spansy(); b++){
			for (int c = 0; c < solid.spansz(); c++){

				//temporary array for dense results of single knot span
				double *values = new double[9*solid.pspantot() * solid.pspantot()];

				Element_BSplineSolid::K_buildSpan(solid,
					D,
					a, ptx, raw_weights_x, dbx, db_1x,
					b, pty, raw_weights_y, dby, db_1y,
					c, ptz, raw_weights_z, dbz, db_1z,
					values);

#pragma omp critical
				{
					Element_BSplineSolid::K_addSpan(solid, a, b, c, values, K);
				}

				//for (int i = 0; i < 20; i++){
				//	std::cout << "i: " << i << ", value = " << values[i] << "\n";
				//}



				delete[] values;

			}
		}
		std::cout << "a: " << a << "\n";
	}

	delete[] raw_abiscas_x;
	delete[] raw_abiscas_y;
	delete[] raw_abiscas_z;
	delete[] raw_weights_x;
	delete[] raw_weights_y;
	delete[] raw_weights_z;
	delete[] dbx;
	delete[] db_1x;
	delete[] weights_x;
	delete[] dby;
	delete[] db_1y;
	delete[] weights_y;
	delete[] dbz;
	delete[] db_1z;
	delete[] weights_z;

	return 0;
}


int Element_BSplineSolid::K_addSpan(BSplineSolid &solid, int span_x, int span_y, int span_z, double values[], SuperMatrix *K){
	//add results to M
	int Kindex1, Kindex2, Kindex3, row;
	int index = 0;

	NCformat *store = (NCformat*)K->Store;
	double* nzval = (double*)store->nzval;


	//for rows, cycle over points in knot span (span_x, span_y, span_z)
	for (int i = 0; i <= solid.p(); i++){
		for (int j = 0; j <= solid.q(); j++){
			for (int k = 0; k <= solid.r(); k++){

				row = 3 * ((i + span_x) * solid.coffsetx() + (j + span_y) * solid.coffsety() + (k + span_z * solid.coffsetz()));

				//for columns, find position of columns in compressed row
				int stepz = 3;
				int stepy = (min(solid.nz(), span_z + k + solid.r() + 1) - max(0, span_z + k - solid.r())) * stepz;
				int stepx = (min(solid.ny(), span_y + j + solid.q() + 1) - max(0, span_y + j - solid.q())) * stepy;

				for (int l = min(span_x, solid.p() - i); l <= min(span_x, solid.p() - i) + solid.p(); l++){
					for (int m = min(span_y, solid.q() - j); m <= min(span_y, solid.q() - j) + solid.q(); m++){
						for (int n = min(span_z, solid.r() - k); n <= min(span_z, solid.r() - k) + solid.r(); n++){

							Kindex1 = store->colptr[row] + l * stepx + m * stepy + n * stepz;
							Kindex2 = store->colptr[row + 1] + l * stepx + m * stepy + n * stepz;
							Kindex3 = store->colptr[row + 2] + l * stepx + m * stepy + n * stepz;

							for (int u = 0; u < 3; u++){
								nzval[Kindex1 + u] += values[index + u];
								nzval[Kindex2 + u] += values[index + 3 + u];
								nzval[Kindex3 + u] += values[index + 6 + u];
							}

							index+=9;
						}
					}
				}
			}
		}
	}

	return 0;
}



