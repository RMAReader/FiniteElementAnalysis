#include "stdafx.h"
#include "MatrixFunctions.h"
#include <omp.h>
#include <amp.h>
#include "iostream"



MatrixFunctions::MatrixFunctions()
{
}


MatrixFunctions::~MatrixFunctions()
{
}

/*
int MatrixFunctions::daxpy(DenseVector &x, DenseVector &y, double da){

	for (int i = 0; i < x.n; i++){
		x.values[i] += da * y.values[i];
	}
	return 0;
}

int MatrixFunctions::dscaladd(DenseVector &x, DenseVector &y, double da){

	for (int i = 0; i < x.n; i++){
		x.values[i] = x.values[i] * da + y.values[i];
	}
	return 0;
}

int MatrixFunctions::dcopy(DenseVector &x, DenseVector &y){

	for (int i = 0; i < x.n; i++){
		x.values[i] = y.values[i];
	}
	return 0;
}

double MatrixFunctions::ddot(DenseVector &x, DenseVector &y){
	double ddot = 0;
	for (int i = 0; i < x.n; i++){
		ddot += x.values[i] * y.values[i];
	}
	return ddot;
}

int MatrixFunctions::add(DenseVector &x, DenseVector &y, DenseVector &xy){
	for (int i = 0; i < x.n; i++){
		xy.values[i] = x.values[i] + y.values[i];
	}
	return 0;
}
int MatrixFunctions::subtract(DenseVector &x, DenseVector &y, DenseVector &xy){
	for (int i = 0; i < x.n; i++){
		xy.values[i] = x.values[i] - y.values[i];
	}
	return 0;
}



int MatrixFunctions::product(SparseCRSMatrix &M, DenseVector &v, DenseVector &Mv){

//omp_set_num_threads(8);
//#pragma omp parallel for
	for (int j = 0; j < M.nrows; j++){
		
		Mv.values[j] = 0;
		for (int i = M.rowptr[j]; i < M.rowptr[j + 1]; i++){

			Mv.values[j] += M.values[i] * v.values[M.colind[i]];

		}

	}
	return 0;
}

int MatrixFunctions::product_amp(SparseCRSMatrix &M, DenseVector &v, DenseVector &Mv){

	using namespace concurrency;

	// Create C++ AMP objects.
	array_view<const double, 1> Mvalues(M.nnzero, M.values);
	array_view<const int, 1> Mcolind(M.nnzero, M.colind);
	array_view<const int, 1> Mrowptr(M.nrows+1, M.rowptr);
	array_view<const double, 1> _v(v.n, v.values);
	array_view<double, 1> _Mv(Mv.n, Mv.values);
	_Mv.discard_data();

	parallel_for_each(
		// Define the compute domain, which is the set of threads that are created.
		_Mv.extent,
		// Define the code to run on each thread on the accelerator.
		[=](index<1> j) restrict(amp)
	
	{

		_Mv[j] = 0;
		for (int i = Mrowptr[j]; i < Mrowptr[j + 1]; i++){

			_Mv[j] += Mvalues[i] * _v[Mcolind[i]];

		}

	});
	return 0;
}


int MatrixFunctions::product(SparseCRSMatrix &M, DenseMatrix &X, DenseMatrix &MX){
	
//	omp_set_num_threads(8);
//#pragma omp parallel for
	for (int j = 0; j < M.nrows; j++){
		for (int k = 0; k < X.rows * X.cols; k += X.rows){
			MX.values[k] = 0;
			for (int i = M.rowptr[j]; i < M.rowptr[j + 1]; i++){

				MX.values[k] += M.values[i] * X.values[k + M.colind[i]];

			}
		}
	}
	return 0;
}



int MatrixFunctions::conjugategradient(SparseCRSMatrix &A,DenseVector &x,DenseVector &b,int maxIterations,double tolerance){

	DenseVector r(x.n);
	DenseVector p(x.n);
	DenseVector Ap(x.n);
	double alpha, beta, rnorm;

	//Ap = A * x
	MatrixFunctions::product(A, x, Ap);
	//r = b - Ap
	MatrixFunctions::subtract(b, Ap, r);
	//p = r
	MatrixFunctions::dcopy(p, r);

	for (int j = 0; j < maxIterations; j++){
		rnorm = MatrixFunctions::ddot(r, r);
		std::cout << "iteration " << j << ", rnorm = " << rnorm << "\n";
		if (rnorm < tolerance) return j;

		//Ap = A * x
		MatrixFunctions::product(A, p, Ap);
		alpha = rnorm / MatrixFunctions::ddot(Ap, p);
		//x = x + alpha * p
		MatrixFunctions::daxpy(x, p, alpha);
		//r = r - alpha * Ap
		MatrixFunctions::daxpy(r, Ap, -alpha);
		beta = MatrixFunctions::ddot(r, r) / rnorm;
		//p = r + beta * p
		MatrixFunctions::dscaladd(p, r, beta);
	}
	return -1;
}


//algorithm 9.1 - using preconditioner M
int MatrixFunctions::leftpreconditionedconjugategradient(SparseCRSMatrix &A, SparseCRSMatrix &Minv,DenseVector &x, DenseVector &b, int maxIterations, double tolerance){

	DenseVector r(x.n);
	DenseVector z(x.n);
	DenseVector p(x.n);
	DenseVector Ap(x.n);
	double alpha, beta, rz;

	//Ax = A * x
	MatrixFunctions::product(A, x, Ap);
	//r = b - Ax
	MatrixFunctions::subtract(b, Ap, r);
	//z = Minv * r
	MatrixFunctions::product(Minv, r, z);
	//p = z
	MatrixFunctions::dcopy(p, z);

	for (int j = 0; j < maxIterations; j++){
		rz = MatrixFunctions::ddot(r, z);
		std::cout << "iteration " << j << ", rnorm = " << rz << "\n";
		if (rz < tolerance) return j;

		//Ap = A * x
		MatrixFunctions::product(A, p, Ap);
		alpha = rz / MatrixFunctions::ddot(Ap, p);
		//x = x + alpha * p
		MatrixFunctions::daxpy(x, p, alpha);
		//r = r - alpha * Ap
		MatrixFunctions::daxpy(r, Ap, -alpha);
		//z = Minv * r
		MatrixFunctions::product(Minv, r, z);
		beta = MatrixFunctions::ddot(r, z) / rz;
		//p = z + beta * p
		MatrixFunctions::dscaladd(p, z, beta);
	}
	return -1;
}

*/