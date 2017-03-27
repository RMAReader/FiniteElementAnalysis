#ifndef _BSPLINE_SOLID_H_
#define _BSPLINE_SOLID_H_


//#ifdef FEADLL_EXPORTS
//#define FEADLL_API __declspec(dllexport) 
//#else
//#define FEADLL_API __declspec(dllimport) 
//#endif



class BSplineSolid
{
public:
	BSplineSolid();
	BSplineSolid &BSplineSolid::operator=(BSplineSolid &source);
	BSplineSolid &BSplineSolid::operator=(BSplineSolid *source);
	BSplineSolid(int p, int q, int r, int nx, int ny, int nz);
	BSplineSolid(double x, double y, double z, int p, int q, int r, int nx, int ny, int nz);
	~BSplineSolid();
	void initialise(int p, int q, int r, int nx, int ny, int nz);


	static void olso_split_spans_x(BSplineSolid *input, BSplineSolid *output, int n);
	static void olso_split_spans_y(BSplineSolid *input, BSplineSolid *output, int n);

	static void BSplineSolid::split_knots(int p, double* input_knot, int lknot, int n, double* output_knot, int output_lknot);
	static int BSplineSolid::number_spans(int p, int lknot);

	double *knotx, *knoty, *knotz;	//knot vectors
	double *cx, *cy, *cz;			//control points

	//getter functions
	/*Number of control points in given dimension*/
	int nx();
	int ny();
	int nz();
	/*Order of bspline in x dimension*/
	int p();
	/*Order of bspline in y dimension*/
	int q();
	/*Order of bspline in z dimension*/
	int r();
	/*Length of knot vector in given dimension*/
	int lknotx();
	int lknoty();
	int lknotz();
	/*Number of knot spans in given dimension*/
	int spansx();
	int spansy();
	int spansz();
	/*Total number of control points*/
	int npoints();
	/*Total number of knot spans*/
	int nspans();
	/*Number of control points affecting single knot span in given dimension*/
	int pspanx();
	int pspany();
	int pspanz();
	/*Total control points affecting single knot span*/
	int pspantot();
	/*Offsets for indexing control point (i, j, k): c_ijk = c[coffsetx*i + coffsety*j + coffsetz * k]*/
	int coffsetx();
	int coffsety();
	int coffsetz();



private:
	int _nx, _ny, _nz;					//size of control net in x, y, z dimensions
	int _p, _q, _r;						//order of bspline in x, y, z dimensions
	int _lknotx, _lknoty, _lknotz;		//knot vector lengths
	int _spansx, _spansy, _spansz;		//number of knot spans in each dimension
	int _npoints, _nspans;		//total number of control pointsand knot spans
	int _pspanx, _pspany, _pspanz;
	int _pspantot;
	int _coffsetx, _coffsety, _coffsetz;
};


#endif _BSPLINE_SOLID_H_