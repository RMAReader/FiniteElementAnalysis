#pragma once



#ifdef FEADLL_EXPORTS
#define FEADLL_API __declspec(dllexport) 
#else
#define FEADLL_API __declspec(dllimport) 
#endif


class FEADLL_API MatrixFunctions
{
public:
	MatrixFunctions();
	~MatrixFunctions();


	/*

	//Level 1 Blas: x = x + da*y
	static int daxpy(DenseVector &x, DenseVector &y, double da);
	static int dscaladd(DenseVector &x, DenseVector &y, double da);
	static int dcopy(DenseVector &x, DenseVector &y);
	static double MatrixFunctions::ddot(DenseVector &x, DenseVector &y);

	static int add(DenseVector &x, DenseVector &y, DenseVector &xy);
	static int subtract(DenseVector &x, DenseVector &y, DenseVector &xy);

	static int product(SparseCRSMatrix &M, DenseVector &v, DenseVector &Mv);
	static int product_amp(SparseCRSMatrix &M, DenseVector &v, DenseVector &Mv);

	static int product(SparseCRSMatrix &M, DenseMatrix &X, DenseMatrix &MX);
	
	//solves Ax = b using the conjugate gradient method
	static int MatrixFunctions::conjugategradient(SparseCRSMatrix &A, DenseVector &x, DenseVector &b, int maxIterations, double tolerance);
	static int MatrixFunctions::leftpreconditionedconjugategradient(SparseCRSMatrix &A, SparseCRSMatrix &Minv, DenseVector &x, DenseVector &b, int maxIterations, double tolerance);
	
	
	*/

};

