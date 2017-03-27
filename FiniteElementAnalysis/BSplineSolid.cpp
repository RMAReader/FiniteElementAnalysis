#include "stdafx.h"
#include "BSplineSolid.h"
#include "math.h"
#include "geometry.h"

BSplineSolid::BSplineSolid()
{
}


/*
Assignment operator
*/
BSplineSolid &BSplineSolid::operator=(BSplineSolid &source)
{
	this->initialise(source.p(), source.q(), source.r(), source.nx(), source.ny(), source.nz());
	for (int i = 0; i < _lknotx; i++){
		knotx[i] = source.knotx[i];
	}
	for (int j = 0; j < _lknoty; j++){
		knoty[j] = source.knoty[j];
	}
	for (int k = 0; k < _lknotz; k++){
		knotz[k] = source.knotz[k];
	}
	int n = source.nx()*source.ny()*source.nz();
	for (int i = 0; i < n; i++){
		cx[i] = source.cx[i];
		cy[i] = source.cy[i];
		cz[i] = source.cz[i];
	}
	return *this;
}


BSplineSolid &BSplineSolid::operator=(BSplineSolid *source)
{
	this->initialise(source->p(), source->q(), source->r(), source->nx(), source->ny(), source->nz());
	for (int i = 0; i < _lknotx; i++){
		knotx[i] = source->knotx[i];
	}
	for (int j = 0; j < _lknoty; j++){
		knoty[j] = source->knoty[j];
	}
	for (int k = 0; k < _lknotz; k++){
		knotz[k] = source->knotz[k];
	}
	int n = source->nx()*source->ny()*source->nz();
	for (int i = 0; i < n; i++){
		cx[i] = source->cx[i];
		cy[i] = source->cy[i];
		cz[i] = source->cz[i];
	}
	return *this;
}

/*
Default contructor creates a cuboid of given size, order and number of points
with control points distributed in a cubic lattice
*/
BSplineSolid::BSplineSolid(double x, double y, double z, int p, int q, int r, int nx, int ny, int nz)
{
	this->initialise(p, q, r, nx, ny, nz);
	
	//default open uniform knot vectors
	for (int i = 0; i < _lknotx; i++){
		knotx[i] = fmin(fmax((double)(i - p) / (nx - p), 0), 1);
	}
	for (int j = 0; j < _lknoty; j++){
		knoty[j] = fmin(fmax((double)(j - q) / (ny - q), 0), 1);
	}
	for (int k = 0; k < _lknotz; k++){
		knotz[k] = fmin(fmax((double)(k - r) / (nz - r), 0), 1);
	}
	int point = 0;
	for (int i = 0; i < nx; i++){
		for (int j = 0; j < ny; j++){
			for (int k = 0; k < nz; k++, point++){
				cx[point] = (double)i * x / (_nx - _p + 1);//(double)i*x / (nx - 1);
				cy[point] = (double)j * y / (_ny - _q + 1);//(double)j*y / (ny - 1);
				cz[point] = (double)k * z / (_nz - _r + 1);//(double)k*z / (nz - 1);
			}
		}
	}
}

BSplineSolid::BSplineSolid(int p, int q, int r, int nx, int ny, int nz)
{
	BSplineSolid::BSplineSolid(1, 1, 1, p, q, r, nx, ny, nz);
}

void BSplineSolid::initialise(int p, int q, int r, int nx, int ny, int nz){
	_p = p;
	_q = q;
	_r = r;
	_nx = nx;
	_ny = ny;
	_nz = nz;

	_lknotx = _nx + _p + 1;
	_lknoty = _ny + _q + 1;
	_lknotz = _nz + _r + 1;

	_spansx = _nx - _p;
	_spansy = _ny - _q;
	_spansz = _nz - _r;

	_npoints = _nx*_ny*_nz;
	_nspans = _spansx*_spansy*_spansz;

	_pspanx = _p + 1;
	_pspany = _q + 1;
	_pspanz = _r + 1;
	_pspantot = _pspanx*_pspany*_pspanz;

	_coffsetx = _ny *_nz;
	_coffsety = _nz;
	_coffsetz = 1;

	knotx = new double[_lknotx];
	knoty = new double[_lknoty];
	knotz = new double[_lknotz];

	int n = nx*ny*nz;
	cx = new double[n];
	cy = new double[n];
	cz = new double[n];
}


BSplineSolid::~BSplineSolid()
{
	delete[] knotx;
	delete[] knoty;
	delete[] knotz;
	delete[] cx;
	delete[] cy;
	delete[] cz;
}

int BSplineSolid::nx(){ return _nx; }
int BSplineSolid::ny(){ return _ny; }
int BSplineSolid::nz(){ return _nz; }
int BSplineSolid::p(){ return _p; }
int BSplineSolid::q(){ return _q; }
int BSplineSolid::r(){ return _r; }
int BSplineSolid::lknotx(){ return _lknotx; }
int BSplineSolid::lknoty(){ return _lknoty; }
int BSplineSolid::lknotz(){ return _lknotz; }
int BSplineSolid::spansx(){ return _spansx; }
int BSplineSolid::spansy(){ return _spansy; }
int BSplineSolid::spansz(){ return _spansz; }
int BSplineSolid::npoints(){ return _npoints; }
int BSplineSolid::nspans(){ return _nspans; }
int BSplineSolid::pspanx(){ return _pspanx; }
int BSplineSolid::pspany(){ return _pspany; }
int BSplineSolid::pspanz(){ return _pspanz; }
int BSplineSolid::pspantot(){ return _pspantot; }
int BSplineSolid::coffsetx(){ return _coffsetx; }
int BSplineSolid::coffsety(){ return _coffsety; }
int BSplineSolid::coffsetz(){ return _coffsetz; }



/*
splits each knot span in x dimension into n spans
*/
void BSplineSolid::olso_split_spans_x(BSplineSolid *input, BSplineSolid *output, int factor){
	if (factor <= 1)return;

	int p = input->p();
	
	int new_nx = input->nx() + input->spansx() * (factor - 1);
	double* cx = new double[input->nx()];
	double* cy = new double[input->nx()];
	double* cz = new double[input->nx()];
	double* dx = new double[new_nx];
	double* dy = new double[new_nx];
	double* dz = new double[new_nx];

	output->initialise(input->p(), input->q(), input->r(), new_nx, input->ny(), input->nz());
	for (int i = 0; i < input->lknoty(); i++){
		output->knoty[i] = input->knoty[i];
	}
	for (int i = 0; i < input->lknotz(); i++){
		output->knotz[i] = input->knotz[i];
	}
	split_knots(input->p(), input->knotx, input->lknotx(), factor, output->knotx, output->lknotx());


	int index;
	for (int j = 0; j < input->ny(); j++){
		for (int k = 0; k < input->nz(); k++){

			for (int i = 0; i < input->nx(); i++){
				index = i*input->coffsetx() + j*input->coffsety() + k*input->coffsetz();
				cx[i] = input->cx[index];
				cy[i] = input->cx[index];
				cz[i] = input->cx[index];
			}
			
			geometry::olso_insertion(input->nx() - 1, cx, input->p() + 1, input->knotx, output->knotx, output->lknotx(), dx);
			geometry::olso_insertion(input->nx() - 1, cy, input->p() + 1, input->knotx, output->knotx, output->lknotx(), dy);
			geometry::olso_insertion(input->nx() - 1, cz, input->p() + 1, input->knotx, output->knotx, output->lknotx(), dz);
			
			for (int i = 0; i < new_nx; i++){
				index = i*input->coffsetx() + j*input->coffsety() + k*input->coffsetz();
				output->cx[index] = dx[i];
				output->cy[index] = dy[i];
				output->cz[index] = dz[i];
			}

		}
	}

	delete[] cx;
	delete[] cy;
	delete[] cz;
	delete[] dx;
	delete[] dy;
	delete[] dz;

}
/*
splits each knot span in y dimension into n spans
*/
void BSplineSolid::olso_split_spans_y(BSplineSolid *input, BSplineSolid *output, int factor){
	if (factor <= 1)return;

	int q = input->q();

	int new_ny = input->ny() + input->spansy() * (factor - 1);
	double* cx = new double[input->ny()];
	double* cy = new double[input->ny()];
	double* cz = new double[input->ny()];
	double* dx = new double[new_ny];
	double* dy = new double[new_ny];
	double* dz = new double[new_ny];

	output->initialise(input->p(), input->q(), input->r(), input->nx(), new_ny, input->nz());
	for (int i = 0; i < input->lknotx(); i++){
		output->knotx[i] = input->knotx[i];
	}
	for (int i = 0; i < input->lknotz(); i++){
		output->knotz[i] = input->knotz[i];
	}
	split_knots(input->q(), input->knoty, input->lknoty(), factor, output->knoty, output->lknoty());


	int index;
	for (int i = 0; i < input->nx(); i++){
		for (int k = 0; k < input->nz(); k++){
			
			for (int j = 0; j < input->ny(); j++){
				index = i*input->coffsetx() + j*input->coffsety() + k*input->coffsetz();
				cx[j] = input->cx[index];
				cy[j] = input->cx[index];
				cz[j] = input->cx[index];
			}

			geometry::olso_insertion(input->ny() - 1, cx, input->q() + 1, input->knoty, output->knoty, output->lknoty(), dx);
			geometry::olso_insertion(input->ny() - 1, cy, input->q() + 1, input->knoty, output->knoty, output->lknoty(), dy);
			geometry::olso_insertion(input->ny() - 1, cz, input->q() + 1, input->knoty, output->knoty, output->lknoty(), dz);

			for (int j = 0; j < new_ny; j++){
				index = i*input->coffsetx() + j*input->coffsety() + k*input->coffsetz();
				output->cx[index] = dx[j];
				output->cy[index] = dy[j];
				output->cz[index] = dz[j];
			}

		}
	}

	delete[] cx;
	delete[] cy;
	delete[] cz;
	delete[] dx;
	delete[] dy;
	delete[] dz;

}



	/*creates a new knot vector with each span split evenly into n spans*/
void BSplineSolid::split_knots(int p, double* input_knot, int lknot, int n, double* output_knot, int output_lknot){
	if (n < 1) return;
	int spans = number_spans(p, lknot);
	if (spans < 1)return;

	/*output_lknot = lknot + spans * (n - 1);
	output_knot = new double[output_lknot];*/

	for (int i = 0; i <= p; i++){
			output_knot[i] = input_knot[i];
			output_knot[output_lknot - i - 1] = input_knot[lknot - i - 1];
	}
		for (int i = 0; i < spans; i++){
			for (int j = 0; j < n; j++){
				output_knot[p + n * i + j] = (double)j / n * input_knot[i + 1] + (1 - (double)j / n)*input_knot[i];
			}
		}

}

	/*Returns the number of knot spans covered by a knot vector*/
int BSplineSolid::number_spans(int p, int lknot){
	return lknot - 2 * p - 1;
}