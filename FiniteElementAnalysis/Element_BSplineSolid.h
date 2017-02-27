#ifndef _ELEMENT_BSPLINE_SOLID_H_
#define _ELEMENT_BSPLINE_SOLID_H_

//#ifdef FEADLL_EXPORTS
//#define FEADLL_API __declspec(dllexport) 
//#else
//#define FEADLL_API __declspec(dllimport) 
//#endif

#include "BSplineSolid.h"
#include "slu_ddefs.h"


class Element_BSplineSolid
{
public:
	Element_BSplineSolid(BSplineSolid &solid, double* D, double density);
	~Element_BSplineSolid();

	
	//defines density of material
	double density;

	//defines elasticity of material. Array contains unique elements of D as {D00, D01, D02, D11, D12, D22, D33, D44, D55}
	double D[9];

	//defines geometry of solid
	BSplineSolid solid;	

	//returns a reference to the inertia matrix in compressed row storage
	SuperMatrix &getM();

	//returns a reference to the stiffness matrix in compressed row storage
	SuperMatrix &getK();

	//returns mass of element
	double get_mass();

public: SuperMatrix M;
public: SuperMatrix K;



	static int buildBasisFunctions(double knot[], int lknot, int p, double raw_abiscas[], double  raw_weights[], int pts, double db[], double db_1[], double weights[]);

	static int getGaussPoints(double knot[], int lknot, int p, int pts, double raw_abiscas[], double  raw_weights[], double abiscas[], double weights[],int lresult);


	static void Element_BSplineSolid::ScaleWeights(double lower, double upper, int n, double raw_weights[], double weights[]);

	static double Element_BSplineSolid::Determinant(double input[]);

	static void Element_BSplineSolid::DetermnantInverse(double M[], double inv[], double &det);

	static void Element_BSplineSolid::Product(double M[], double v[], double Mv[]);


	static int Element_BSplineSolid::buildJacobianDet(BSplineSolid &solid,
		int span_x, int ptx, double db_x[], double db_x1[],
		int span_y, int pty, double db_y[], double db_y1[],
		int span_z, int ptz, double db_z[], double db_z1[],
		double JacobianDet[]);
	
	static void Element_BSplineSolid::initialise_SparseCRSMatrix(BSplineSolid &solid, SuperMatrix*);

	static int Element_BSplineSolid::M_build(BSplineSolid &solid, double density, SuperMatrix*);

	static int Element_BSplineSolid::M_addSpan(BSplineSolid &solid, int span_x, int span_y, int span_z, double values[], SuperMatrix* );

	static int Element_BSplineSolid::M_buildSpan(BSplineSolid &solid, double density,
		int span_x, int ptx, double db_x[], double db_x1[], double weights_x[],
		int span_y, int pty, double db_y[], double db_y1[], double weights_y[],
		int span_z, int ptz, double db_z[], double db_z1[], double weights_z[],
		double values[]);


	static int Element_BSplineSolid::K_build(BSplineSolid &solid, double* D, SuperMatrix*);

	static int Element_BSplineSolid::K_buildSpan(BSplineSolid &solid, double* D,
		int span_x, int ptx, double db_x[], double db_x1[], double weights_x[],
		int span_y, int pty, double db_y[], double db_y1[], double weights_y[],
		int span_z, int ptz, double db_z[], double db_z1[], double weights_z[],
		double values[]);

	static int Element_BSplineSolid::K_addSpan(BSplineSolid &solid, int span_x, int span_y, int span_z, double values[], SuperMatrix*);




};




#endif _ELEMENT_BSPLINE_SOLID_H_