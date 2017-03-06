// TestFiniteElementAnalysis.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "iostream"
#include <math.h>
//#include "HarwellBoeing.h"
#include "BSplineSolid.h"
#include "Element_BSplineSolid.h"
#include "gauss_legendre.h"
#include "MatrixFunctions.h"
#include <chrono>
#include "arpack_shiftinvert.h"
#include "slu_ddefs.h"
#include "FileReader.h"
#include "Optimizer.h"
#include "bspline_utils.h"
#include "violin_model.h"
#include "toolpath_base.h"

//void TestHarwellBoeing();
bool TestDeBoor(bool verbose);
bool TestBasisFunctions(bool verbose);
void TestM();
void TestK();
void TestDeterminants();
void TestMatrixFunctions();
void TestSuperLU();
void TestEigenSolver();
//bool TestFileLoader(bool verbose);
bool TestFileLoader_LoadModelJava(bool verbose);
bool TestFileLoader_LoadViolinXML(bool verbose);
void TestOptimizer();
bool TestOlsoAlgorithm(bool verbose);
bool TestSurface(bool verbose);
bool TestCurve_InsertKnots(bool verbose);
bool TestCurve_LineIntersection(bool verbose);
bool TestCurve_CurveIntersection(bool verbose);
bool TestToolPath3(bool verbose);
bool TestTriangle(bool verbose);

bool AreEqual(double p, double q, double error){
	if (abs(p - q) > error){ return false; }
	return true;
}
bool NotEqual(double p, double q, double error){
	if (abs(p - q) <= error){ return false; }
	return true;
}

int _tmain(int argc, _TCHAR* argv[])
{
	int TestsPassed = 0;
	int TestsFailed = 0;
	bool verbose = false;


	if (TestDeBoor(verbose))
	{
		cout << "TestDeBoor passed" << endl; TestsPassed++;
	}
	else{
		cout << "TestDeBoor failed" << endl; TestsFailed++;
	}

	if (TestBasisFunctions(verbose))
	{
		cout << "TestBasisFunctions passed" << endl; TestsPassed++;
	}
	else{
		cout << "TestBasisFunctions failed" << endl; TestsFailed++;
	}

	//TestHarwellBoeing();
	//TestM();
	//TestK();
	//TestDeterminants();
	//TestMatrixFunctions();
	//TestSuperLU();
	//TestEigenSolver();
	//TestFileLoader();
	//TestFileLoader2();
	//TestOptimizer();
	
	if (TestOlsoAlgorithm(verbose))
	{
		cout << "TestOlsoAlgorithm passed" << endl; TestsPassed++;
	}
	else{
		cout << "TestOlsoAlgorithm failed" << endl; TestsFailed++;
	}


	if (TestSurface(verbose))
	{
		cout << "TestSurface passed" << endl; TestsPassed++;
	}
	else{
		cout << "TestSurface failed" << endl; TestsFailed++;
	}

	if (TestCurve_InsertKnots(verbose))
	{
		cout << "TestCurve_InsertKnots passed" << endl; TestsPassed++;
	}
	else{
		cout << "TestCurve_InsertKnots failed" << endl; TestsFailed++;
	}

	if (TestCurve_LineIntersection(verbose))
	{
		cout << "TestCurve_Intersection passed" << endl; TestsPassed++;
	}
	else{ 
		cout << "TestCurve_Intersection failed" << endl; TestsFailed++; 
	}

	if (TestCurve_CurveIntersection(verbose))
	{
		cout << "TestCurve_CurveInterction passed" << endl; TestsPassed++;
	}
	else{
		cout << "TestCurve_CurveInterction failed" << endl; TestsFailed++;
	}

	if (TestTriangle(verbose))
	{
		cout << "TestTriangle passed" << endl; TestsPassed++;
	}
	else{
		cout << "TestTriangle failed" << endl; TestsFailed++;
	}


	if (TestFileLoader_LoadModelJava(verbose))
	{
		cout << "TestFileLoader_LoadModelJava passed" << endl; TestsPassed++;
	}
	else{
		cout << "TestFileLoader_LoadModelJava failed" << endl; TestsFailed++;
	}


	if (TestFileLoader_LoadViolinXML(verbose))
	{
		cout << "TestFileLoader_LoadViolinXML passed" << endl; TestsPassed++;
	}
	else{
		cout << "TestFileLoader_LoadViolinXML failed" << endl; TestsFailed++;
	}

	//if (TestToolPath1(verbose))
	//{
	//	cout << "TestToolPath1 passed" << endl; TestsPassed++;
	//}
	//else{
	//	cout << "TestToolPath1 failed" << endl; TestsFailed++;
	//}

	//if (TestToolPath2(verbose))
	//{
	//	cout << "TestToolPath2 passed" << endl; TestsPassed++;
	//}
	//else{
	//	cout << "TestToolPath2 failed" << endl; TestsFailed++;
	//}

	if (TestToolPath3(verbose))
	{
		cout << "TestToolPath3 passed" << endl; TestsPassed++;
	}
	else{
		cout << "TestToolPath3 failed" << endl; TestsFailed++;
	}

	cout << endl;
	cout << "Summary:" << endl;
	cout << "Number of tests passed: " << TestsPassed << endl;
	cout << "Number of tests failed: " << TestsFailed << endl;
	cout << endl;

	system("pause");

	return 0;
}


bool TestDeBoor(bool verbose){
	int i, p, lknot;
	double x;
	double* knot;
	double tol = 0.000001;

	//case 1
	int p1 = 2;
	int lknot1 = 6;
	double knot1[6] = { 0.01, 0.04, 0.09, 0.16, 0.25, 0.36 };
	double evaluation_points[6] = { 0.01, 0.03, 0.04, 0.07, 0.13 ,0.16};
	double expected_values[6] = { 0, 0.166667, 0.375, 0.75, 0.107143, 0.0 };

	if (verbose){
		std::cout << "i=0" << endl;
		std::cout << "p=" << p1 << endl;
		for (int i = 0; i < 5; i++){
			std::cout << "t = " << evaluation_points[i];
			std::cout << "b(t) = " << bspline::deboor_value(0, p1, knot1, lknot1, evaluation_points[i]);
			std::cout << endl;
		}
	}

	for (int i = 0; i < 6; i++){
		x = bspline::deboor_value(0, p1, knot1, lknot1, evaluation_points[i]);
		if (AreEqual(x, expected_values[i], tol) == false){ return false; }
	}



	//std::cout << "Enter i: ";
	//std::cin >> i;

	//std::cout << "Enter p: ";
	//std::cin >> p;

	//std::cout << "Enter length of knot vector: ";
	//std::cin >> lknot;

	//knot = new double[lknot];

	//for (int j = 0; j < lknot; j++){
	//	std::cout << "knot[" << j << "] = ";
	//	std::cin >> knot[j];
	//}

	//std::cout << "Enter x: ";
	//std::cin >> x;

	//std::cout << "deBoor value = " << bspline::deboor_value(i, p, knot, lknot, x) << "\n";

	//for (int j = 0; j < lknot * 10; j++){
	//	x = j * 0.1;
	//	double* values = new double[p];
	//	bspline::deboor_values(i, p, knot, lknot, x, values);
	//	std::cout << "deBoor(";
	//	std::cout << x;
	//	std::cout << ") = ";
	//	std::cout << bspline::deboor_value(i, p, knot, lknot, x);
	//	std::cout << ", ";
	//	std::cout << values[0];
	//	std::cout << ", ";
	//	std::cout << bspline::deboor_derivative(i, p, knot, lknot, x);
	//	std::cout << ", ";
	//	std::cout << values[1];
	//	std::cout << "\n";
	//}


	//if (true){
	//	return true;
	//}
	//else{
	//	return false;
	//}
	

}

//void TestHarwellBoeing(){
//	std::string filename = "TestMatrix.rsa";
//	int nrow = 5;
//	int ncol = 5;
//	int nnzero = 9;
//	int rowind[] = { 1, 1, 2, 2, 3, 3, 4, 4, 5 };
//	int colptr[] = { 1, 3, 5, 7, 9 };
//	double values[] = { 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9 };
//
//	HarwellBoeing::Save(filename, nrow, ncol, nnzero, values, rowind, colptr);
//
//}
//
//
//
//
bool TestBasisFunctions(bool verbose){

	double tol = 0.0000001;

	//8-point Gauss points in [-1,1]
	double raw_abiscas[8];// = { -0.1834346424956498049394761, -0.5255324099163289858177390, -0.7966664774136267395915539, -0.9602898564975362316835609, 0.1834346424956498049394761, 0.5255324099163289858177390, 0.7966664774136267395915539, 0.9602898564975362316835609 };
	double raw_weights[8];// = { 0.3626837833783619829651504, 0.3137066458778872873379622, 0.2223810344533744705443560, 0.1012285362903762591525314, 0.3626837833783619829651504, 0.3137066458778872873379622, 0.2223810344533744705443560, 0.1012285362903762591525314 };
	int pts = 8;

	gauss_legendre_tbl(pts, raw_abiscas, raw_weights, 1E-15);

	double expected_raw_abiscas[8] = { 0.1834346424956498049394761, 0.5255324099163289858177390, 0.7966664774136267395915539, 0.9602898564975362316835609, -0.1834346424956498049394761, -0.5255324099163289858177390, -0.7966664774136267395915539, -0.9602898564975362316835609 };
	double expected_raw_weights[8] = { 0.3626837833783619829651504, 0.3137066458778872873379622, 0.2223810344533744705443560, 0.1012285362903762591525314, 0.3626837833783619829651504, 0.3137066458778872873379622, 0.2223810344533744705443560, 0.1012285362903762591525314 };

	for (int i = 0; i < 8; i++){
		if (AreEqual(raw_abiscas[i], expected_raw_abiscas[i], tol) == false){ return false; }
		if (AreEqual(raw_weights[i], expected_raw_weights[i], tol) == false){ return false; }
	}

	double knot[] = { 0, 1, 2, 3, 4, 5, 6, 7 };
	int lknot = 8;
	int p = 2;

	int nSpans = (lknot - 2 * p-1);			//number of knot spans
	int nCPoints = nSpans * (p + 1);		//number of spans * control point which affect each span
	int nIPoints = nCPoints * pts;			//number of spans * control point which affect each span * integration points per span

	double* db = new double[nIPoints];		//ouput array
	double* db_1 = new double[nIPoints];	//ouput array
	double* weights = new double[nIPoints];	//ouput array

	int res1 = Element_BSplineSolid::buildBasisFunctions(knot, lknot, p, raw_abiscas, raw_weights, pts, db, db_1, weights);

	if (verbose){
		for (int i = 0; i < nIPoints; i++){
			std::cout << i << ": " << db[i] << ", " << db_1[i] << ", " << weights[i] << "\n";
		}
	}
	return true;
}

void TestDeterminants(){

	double M[9] = { 1, 2, 3,
		4, 6, 3,
		6, 9, -6 };

	double det1 = Element_BSplineSolid::Determinant(M);

	double Minv[9];
	double det2;

	Element_BSplineSolid::DetermnantInverse(M,Minv,det2);

	std::cout <<  "M:\n";
	for (int i = 0; i < 9; i++){
				std::cout << i << ": " << M[i] << "\n";
	}
	std::cout << "Minv:\n";
	for (int i = 0; i < 9; i++){
		std::cout << i << ": " << Minv[i] << "\n";
	}
	std::cout << "det1 = " << det1 << "\n";
	std::cout << "det2 = " << det2 << "\n";

}



//calculate M for a 2-order 2x2x2 cube, i.e. 4x4x4 control points
void TestM(){
	std::cout << "Beginning M test\n";

	double density = 2459.9;
	BSplineSolid solid(0.002333, 0.002889, 0.003914, 2, 2, 2, 50, 50, 50);

	std::cout << "p = " << solid.p() << "\n";
	std::cout << "q = " << solid.q() << "\n";
	std::cout << "r = " << solid.r() << "\n";
	std::cout << "nx = " << solid.nx() << "\n";
	std::cout << "ny = " << solid.ny() << "\n";
	std::cout << "nz = " << solid.nz() << "\n";

	std::cout << "M dimension = " << 3 * solid.npoints() << "\n";
	

	SuperMatrix M;

	Element_BSplineSolid::M_build(solid, density, &M);
		
	std::cout << "Test M compete\n";
	NCformat *store = (NCformat*) M.Store;
	double* nzval = (double*)store->nzval;

	std::cout << "nnzero = " << store->nnz<< "\n";
	std::cout << "density = " << (double)store->nnz / M.nrow / M.nrow << "\n";
	std::cout << "memory = " << (store->nnz * 12 + M.nrow * 4) / 1000000 << "MB\n";

	int row = 0;
	for (int i = 0; i < 10; i++){
		if (i == store->colptr[row+1]){ row++; }
		std::cout << "Row: " << row << ", ";
		std::cout << "Column: " << store->rowind[i] << ", ";
		std::cout << "Value: " << nzval[i] << "\n";
	}
	for (int i = store->nnz / 2 - 5; i < store->nnz / 2 + 5; i++){
		if (i == store->colptr[row + 1]){ row++; }
		std::cout << "Row: " << row << ", ";
		std::cout << "Column: " << store->rowind[i] << ", ";
		std::cout << "Value: " << nzval[i] << "\n";
	}
	for (int i = store->nnz-10; i < store->nnz; i++){
		if (i == store->colptr[row + 1]){ row++; }
		std::cout << "Row: " << row << ", ";
		std::cout << "Column: " << store->rowind[i] << ", ";
		std::cout << "Value: " << nzval[i] << "\n";
	}


}


void TestK(){
	std::cout << "Beginning K test\n";

	//glass
	double density = 2459.9;
	//double D1[6][6] = { { 82.0407E9, 23.5666E9, 23.5666E9,         0,         0,         0 },
	//					{ 23.5666E9, 82.0407E9, 23.5666E9,         0,         0,         0 },
	//					{ 23.5666E9, 23.5666E9, 82.0407E9,         0,         0,         0 },
	//					{         0,         0,         0, 29.2371E9,         0,         0 },
	//					{         0,         0,         0,         0, 29.2371E9,         0 },
	//					{         0,         0,         0,         0,         0, 29.2371E9 }};

	double D1[9] = { 82.0407E9, 23.5666E9, 23.5666E9, 82.0407E9, 23.5666E9, 82.0407E9, 29.2371E9, 29.2371E9, 29.2371E9 };

	BSplineSolid solid(0.002333, 0.002889, 0.003914, 2, 2, 2, 30, 30, 30);

	std::cout << "p = " << solid.p() << "\n";
	std::cout << "q = " << solid.q() << "\n";
	std::cout << "r = " << solid.r() << "\n";
	std::cout << "nx = " << solid.nx() << "\n";
	std::cout << "ny = " << solid.ny() << "\n";
	std::cout << "nz = " << solid.nz() << "\n";

	std::cout << "dimension = " << 3 * solid.npoints() << "\n";

	Element_BSplineSolid element(solid,D1,density);

	SuperMatrix K = element.getK();
	SuperMatrix M = element.getM();

	std::cout << "Test K compete\n";
	NCformat *store = (NCformat*)K.Store;
	double* nzval = (double*)store->nzval;

	std::cout << "nnzero = " << store->nnz << "\n";
	std::cout << "density = " << (double)store->nnz / K.nrow / K.nrow << "\n";
	std::cout << "memory = " << (store->nnz * 12 + K.nrow * 4) / 1000000 << "MB\n";

	int row = 0;
	for (int i = 0; i < 10; i++){
		if (i == store->colptr[row + 1]){ row++; }
		std::cout << "Row: " << row << ", ";
		std::cout << "Column: " << store->rowind[i] << ", ";
		std::cout << "Value: " << nzval[i] << "\n";
	}
	for (int i = store->nnz / 2 - 5; i < store->nnz / 2 + 5; i++){
		if (i == store->colptr[row + 1]){ row++; }
		std::cout << "Row: " << row << ", ";
		std::cout << "Column: " << store->rowind[i] << ", ";
		std::cout << "Value: " << nzval[i] << "\n";
	}
	for (int i = store->nnz - 10; i < store->nnz; i++){
		if (i == store->colptr[row + 1]){ row++; }
		std::cout << "Row: " << row << ", ";
		std::cout << "Column: " << store->rowind[i] << ", ";
		std::cout << "Value: " << nzval[i] << "\n";
	}


}


/*
void TestMatrixFunctions(){

	//glass
	double density = 2459.9;
	double D1[6][6] = { { 82.0407E9, 23.5666E9, 23.5666E9, 0, 0, 0 },
	{ 23.5666E9, 82.0407E9, 23.5666E9, 0, 0, 0 },
	{ 23.5666E9, 23.5666E9, 82.0407E9, 0, 0, 0 },
	{ 0, 0, 0, 29.2371E9, 0, 0 },
	{ 0, 0, 0, 0, 29.2371E9, 0 },
	{ 0, 0, 0, 0, 0, 29.2371E9 } };


	BSplineSolid solid(0.002333, 0.002889, 0.003914, 2, 2, 2, 10, 10, 10);

	std::cout << "p = " << solid.p() << "\n";
	std::cout << "q = " << solid.q() << "\n";
	std::cout << "r = " << solid.r() << "\n";
	std::cout << "nx = " << solid.nx() << "\n";
	std::cout << "ny = " << solid.ny() << "\n";
	std::cout << "nz = " << solid.nz() << "\n";

	std::cout << "dimension = " << 3 * solid.npoints() << "\n";

	Element_BSplineSolid element(solid, D1, density);

	SparseCRSMatrix K = element.getK();
	SparseCRSMatrix M = element.getM();
	std::cout << "K.nnzero = " << K.nnzero << "\n";
	std::cout << "K.nrows = " << K.nrows << "\n";
	std::cout << "K.ncols = " << K.ncols << "\n";



	DenseVector v(K.nrows, 1);

	DenseVector Kv(K.nrows, 0);

	std::cout << "Begin matrix multiply: OpenMP\n";
	auto t1 = std::chrono::high_resolution_clock::now();
	MatrixFunctions::product(K, v, Kv);
	auto t2 = std::chrono::high_resolution_clock::now();
	
	std::cout << "Completed matrix multiply in ";
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
	std::cout << " milliseconds\n";

	//std::cout << "Begin matrix multiply: amp.h\n";
	//t1 = std::chrono::high_resolution_clock::now();
	//MatrixFunctions::product_amp(K, v, Kv);
	//t2 = std::chrono::high_resolution_clock::now();

	//std::cout << "Completed matrix multiply in ";
	//std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
	//std::cout << " milliseconds\n";


	DenseMatrix X(K.nrows, 20, 1);
	DenseMatrix KX(K.nrows, 20, 0);
	std::cout << "Begin matrix multiply: OpenMP\n";
	t1 = std::chrono::high_resolution_clock::now();
	MatrixFunctions::product(K, X, KX);
	t2 = std::chrono::high_resolution_clock::now();

	std::cout << "Completed matrix multiply in ";
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
	std::cout << " milliseconds\n";


	DenseVector b(K.nrows, rand());
	DenseVector x(K.nrows, 0);
	DenseVector y(K.nrows, 0);
	int maxIterations = 500;
	double tolerance = 0.001;
	int result;


	std::cout << "Begin conjugate gradient solve M\n";
	t1 = std::chrono::high_resolution_clock::now();
	result = MatrixFunctions::conjugategradient(M, y, b, maxIterations, tolerance);
	t2 = std::chrono::high_resolution_clock::now();
	std::cout << "Conjugate gradient result = " << result << "\n";
	std::cout << "Completed conjugate gradient solve M in ";
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
	std::cout << " milliseconds\n";


	double shift = 39000000;
	for (int i = 0; i < K.nnzero; i++){
		K.values[i] -= shift * M.values[i];
	}


	std::cout << "Begin conjugate gradient solve (K-aM)\n";
	t1 = std::chrono::high_resolution_clock::now();
	result = MatrixFunctions::conjugategradient(K, x, b, maxIterations, tolerance);
	t2 = std::chrono::high_resolution_clock::now();
	std::cout << "Conjugate gradient result = " << result << "\n";
	std::cout << "Completed conjugate gradient solve (K-aM) in ";
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
	std::cout << " milliseconds\n";

	//diagonal condictioner
	SparseCRSMatrix Conditioner;
	Conditioner.ncols = K.ncols;
	Conditioner.nrows = K.nrows;
	Conditioner.nnzero = K.nrows;
	Conditioner.colind = new int[K.ncols];
	Conditioner.rowptr = new int[K.ncols+1];
	Conditioner.values = new double[K.ncols];
	for (int i = 0; i < K.ncols; i++){
		Conditioner.colind[i] = i;
		Conditioner.rowptr[i] = i;
		for (int j = K.rowptr[i]; j < K.rowptr[i + 1]; j++){
			if (K.colind[j] == i)
			{
				Conditioner.values[i] = 1 / K.values[j];
			}
		}
	}
	Conditioner.rowptr[K.ncols] = K.ncols;
	for (int i = 0; i < x.n; i++){
		x.values[i] = 0;
	}

	std::cout << "Begin preconditioned conjugate gradient solve (K-aM)\n";
	t1 = std::chrono::high_resolution_clock::now();
	result = MatrixFunctions::leftpreconditionedconjugategradient(K, Conditioner,x, b, maxIterations, tolerance);
	t2 = std::chrono::high_resolution_clock::now();
	std::cout << "Conjugate gradient result = " << result << "\n";
	std::cout << "Completed conjugate gradient solve (K-aM) in ";
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
	std::cout << " milliseconds\n";




}

*/

void TestSuperLU(){

	std::cout << "Beginning K test\n";

	//glass
	double density = 2459.9;
	//double D1[6][6] = { { 82.0407E9, 23.5666E9, 23.5666E9, 0, 0, 0 },
	//{ 23.5666E9, 82.0407E9, 23.5666E9, 0, 0, 0 },
	//{ 23.5666E9, 23.5666E9, 82.0407E9, 0, 0, 0 },
	//{ 0, 0, 0, 29.2371E9, 0, 0 },
	//{ 0, 0, 0, 0, 29.2371E9, 0 },
	//{ 0, 0, 0, 0, 0, 29.2371E9 } };
	double D1[9] = { 82.0407E9, 23.5666E9, 23.5666E9, 82.0407E9, 23.5666E9, 82.0407E9, 29.2371E9, 29.2371E9, 29.2371E9 };


	BSplineSolid solid(0.002333, 0.002889, 0.003914, 2, 2, 2, 20, 20, 20);

	std::cout << "p = " << solid.p() << "\n";
	std::cout << "q = " << solid.q() << "\n";
	std::cout << "r = " << solid.r() << "\n";
	std::cout << "nx = " << solid.nx() << "\n";
	std::cout << "ny = " << solid.ny() << "\n";
	std::cout << "nz = " << solid.nz() << "\n";

	std::cout << "dimension = " << 3 * solid.npoints() << "\n";

	Element_BSplineSolid element(solid, D1, density);

	SuperMatrix K = element.getK();
	SuperMatrix M = element.getM();

	SuperMatrix A;
	NCformat *Astore;
	double   *a;
	int      *asub, *xa;
	int      *perm_c; /* column permutation vector */
	int      *perm_r; /* row permutations from partial pivoting */
	SuperMatrix L;      /* factor L */
	SCformat *Lstore;
	SuperMatrix U;      /* factor U */
	NCformat *Ustore;
	SuperMatrix B;
	int      nrhs, ldx, info, m, n, nnz;
	double   *xact, *rhs;
	mem_usage_t   mem_usage;
	superlu_options_t options;
	SuperLUStat_t stat;
	FILE      *fp = stdin;

#if ( DEBUGlevel>=1 )
	CHECK_MALLOC("Enter main()");
#endif

	/* Set the default input options:
	options.Fact = DOFACT;
	options.Equil = YES;
	options.ColPerm = COLAMD;
	options.DiagPivotThresh = 1.0;
	options.Trans = NOTRANS;
	options.IterRefine = NOREFINE;
	options.SymmetricMode = NO;
	options.PivotGrowth = NO;
	options.ConditionNumber = NO;
	options.PrintStat = YES;
	*/
	set_default_options(&options);

	/* Read the matrix in Harwell-Boeing format. */
	//dreadhb(fp, &m, &n, &nnz, &a, &asub, &xa);
	A = K;
	m = K.nrow;
	n = K.ncol;
	NCformat *store = (NCformat*) K.Store;
	nnz = store->nnz;

	//dCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);
	Astore = (NCformat*)A.Store;
	printf("Dimension %dx%d; # nonzeros %d\n", A.nrow, A.ncol, Astore->nnz);

	nrhs = 1;
	if (!(rhs = doubleMalloc(m * nrhs))) ABORT("Malloc fails for rhs[].");
	dCreate_Dense_Matrix(&B, m, nrhs, rhs, m, SLU_DN, SLU_D, SLU_GE);
	xact = doubleMalloc(n * nrhs);
	ldx = n;
	dGenXtrue(n, nrhs, xact, ldx);
	dFillRHS(options.Trans, nrhs, xact, ldx, &A, &B);

	if (!(perm_c = intMalloc(n))) ABORT("Malloc fails for perm_c[].");
	if (!(perm_r = intMalloc(m))) ABORT("Malloc fails for perm_r[].");

	/* Initialize the statistics variables. */
	StatInit(&stat);

	dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);

	if (info == 0) {

		/* This is how you could access the solution matrix. */
		double *sol = (double*)((DNformat*)B.Store)->nzval;

		/* Compute the infinity norm of the error. */
		dinf_norm_error(nrhs, &B, xact);

		Lstore = (SCformat *)L.Store;
		Ustore = (NCformat *)U.Store;
		printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
		printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
		printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz - n);
		printf("FILL ratio = %.1f\n", (float)(Lstore->nnz + Ustore->nnz - n) / nnz);

		dQuerySpace(&L, &U, &mem_usage);
		printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
			mem_usage.for_lu / 1e6, mem_usage.total_needed / 1e6);

	}
	else {
		printf("dgssv() error returns INFO= %d\n", info);
		if (info <= n) { /* factorization completes */
			dQuerySpace(&L, &U, &mem_usage);
			printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
				mem_usage.for_lu / 1e6, mem_usage.total_needed / 1e6);
		}
	}

	if (options.PrintStat) StatPrint(&stat);
	StatFree(&stat);

	SUPERLU_FREE(rhs);
	SUPERLU_FREE(xact);
	SUPERLU_FREE(perm_r);
	SUPERLU_FREE(perm_c);
	Destroy_CompCol_Matrix(&A);
	Destroy_SuperMatrix_Store(&B);
	Destroy_SuperNode_Matrix(&L);
	Destroy_CompCol_Matrix(&U);

#if ( DEBUGlevel>=1 )
	CHECK_MALLOC("Exit main()");
#endif
	

}



void TestEigenSolver(){

	std::cout << "Beginning Eigensolver test\n";

	//glass
	double density = 2459.9;
	//double D1[6][6] = { { 82.0407E9, 23.5666E9, 23.5666E9, 0, 0, 0 },
	//{ 23.5666E9, 82.0407E9, 23.5666E9, 0, 0, 0 },
	//{ 23.5666E9, 23.5666E9, 82.0407E9, 0, 0, 0 },
	//{ 0, 0, 0, 29.2371E9, 0, 0 },
	//{ 0, 0, 0, 0, 29.2371E9, 0 },
	//{ 0, 0, 0, 0, 0, 29.2371E9 } };

	double D1[9] = { 82.0407E9, 23.5666E9, 23.5666E9, 82.0407E9, 23.5666E9,  82.0407E9,  29.2371E9, 29.2371E9,  29.2371E9  };


	BSplineSolid solid(0.002333, 0.002889, 0.003914, 2, 2, 2, 20, 20, 20);

	std::cout << "p = " << solid.p() << "\n";
	std::cout << "q = " << solid.q() << "\n";
	std::cout << "r = " << solid.r() << "\n";
	std::cout << "nx = " << solid.nx() << "\n";
	std::cout << "ny = " << solid.ny() << "\n";
	std::cout << "nz = " << solid.nz() << "\n";

	std::cout << "dimension = " << 3 * solid.npoints() << "\n";

	Element_BSplineSolid element(solid, D1, density);

	arpack_shiftinvert solver(element.getM(), element.getK());
	solver.shift(0);
	solver.factorizeLU();
	solver.solve();


	density *= 1.1;
	Element_BSplineSolid element2(solid, D1, density);

	solver.setM(element2.getM());
	solver.setA(element2.getK());
	solver.solve();
}


bool TestFileLoader(bool verbose){

	if (verbose){ std::cout << "Beginning Loader test\n"; }
	
	char filepath[] = "C:\\Private\\Harris_tuning_dx0";

	double density;
	double D1[9];
	double D[36];
	BSplineSolid solid;
	IO::ReadInput(filepath, &solid, D, density);

	D1[0] = D[0];
	D1[1] = D[1];
	D1[2] = D[2];
	D1[3] = D[8];
	D1[4] = D[9];
	D1[5] = D[14];
	D1[6] = D[21];
	D1[7] = D[28];
	D1[8] = D[35];
	
	std::cout << "p = " << solid.p() << "\n";
	std::cout << "q = " << solid.q() << "\n";
	std::cout << "r = " << solid.r() << "\n";
	std::cout << "nx = " << solid.nx() << "\n";
	std::cout << "ny = " << solid.ny() << "\n";
	std::cout << "nz = " << solid.nz() << "\n";

	std::cout << "dimension = " << 3 * solid.npoints() << "\n";

	Element_BSplineSolid element(solid, D1, 1.1*density);

	arpack_shiftinvert solver(element.getM(), element.getK());
	solver.shift(0);
	solver.factorizeLU();
	solver.solve();


	density /= 1.1;
	Element_BSplineSolid element2(solid, D1, density);

	solver.setM(element2.getM());
	solver.setA(element2.getK());
	solver.factorizeLU();

	solver.solve();

	return true;
}


bool TestFileLoader_LoadModelJava(bool verbose){



	if (verbose){ std::cout << "Beginning Loader test\n"; }

	char filepath[] = "C:\\Users\\Lizzie\\Documents\\GitHub\\FiniteElementAnalysis\\Data\\Perlman Strad belly full v24 high wide v5.dat";

	std::vector<CURVE2F> curves;
	std::vector<SURFACE3F> surfaces;
	bool successful = IO::LoadModelJava(filepath, &curves, &surfaces, verbose);
	
	violin_model* violin = new violin_model();
	violin->ribs = new violin_ribs();
	std::string name = "rib_internal_lower_bout";
	violin->ribs->curves.insert(std::make_pair(name, &(curves[0])));

	char filepath2[] = "C:\\Users\\Lizzie\\Documents\\GitHub\\FiniteElementAnalysis\\Data\\Perlman Strad Violin Model saved.xml";
	IO::SaveModelXML(filepath2, violin, verbose);


	if (successful == true && curves.size() == 16){ 
		return true; 
	}
	else{ return false; }

}



bool TestFileLoader_LoadViolinXML(bool verbose){

	if (verbose){ std::cout << "Beginning Loader test\n"; }

	char filepath[] = "C:\\Users\\Lizzie\\Documents\\GitHub\\Finite-Element-Engine\\FiniteElementAnalysis\\Data\\Perlman Strad Violin Model.xml";
	char filepath2[] = "C:\\Users\\Lizzie\\Documents\\GitHub\\Finite-Element-Engine\\FiniteElementAnalysis\\Data\\Perlman Strad Violin Model saved.xml";

	app_model* model = IO::LoadModelXML(filepath, verbose);
	if (model == nullptr){ return false; }

	IO::SaveModelXML(filepath2, model->violin, verbose);

	return true;
}


void TestOptimizer(){
	std::cout << "Beginning Optimizer test\n";
	////glass
	//double density = 2459.9;
	//double D[6][6] = { { 82.0407E9, 23.5666E9, 23.5666E9, 0, 0, 0 },
	//{ 23.5666E9, 82.0407E9, 23.5666E9, 0, 0, 0 },
	//{ 23.5666E9, 23.5666E9, 82.0407E9, 0, 0, 0 },
	//{ 0, 0, 0, 29.2371E9, 0, 0 },
	//{ 0, 0, 0, 0, 29.2371E9, 0 },
	//{ 0, 0, 0, 0, 0, 29.2371E9 } };

	//BSplineSolid solid(0.002333, 0.002889, 0.003914, 2, 2, 2, 5, 5, 5);

	char filepath[] = "C:\\Private\\Harris_tuning_dx0";
	double density;
	double D1[9];
	double D[36];
	BSplineSolid solid;
	IO::ReadInput(filepath, &solid, D, density);

	std::cout << "p = " << solid.p() << "\n";
	std::cout << "q = " << solid.q() << "\n";
	std::cout << "r = " << solid.r() << "\n";
	std::cout << "nx = " << solid.nx() << "\n";
	std::cout << "ny = " << solid.ny() << "\n";
	std::cout << "nz = " << solid.nz() << "\n";

	//D1[0] = D[0];
	//D1[1] = D[1];
	//D1[2] = D[2];
	//D1[3] = D[7];
	//D1[4] = D[8];
	//D1[5] = D[14];
	//D1[6] = D[21];
	//D1[7] = D[28];
	//D1[8] = D[35];

	//int number_of_modes = 3;
	//int measured_modes[3] = { 1, 2, 5 };
	//double measured_frequencies[3] = { 90.0, 180.0, 360.0 };

	Optimizer optimizer;
	//optimizer.calibrate(&solid, density, D1, measured_modes, measured_frequencies,number_of_modes);



	D1[0] = 10692000000000.00;
	D1[1] = 2165520000.00;
	D1[2] = 165836000.00;
	D1[3] = 949212000.00;
	D1[4] = 53823400.00;
	D1[5] = 320534000.00;
	D1[6] = 5323230.00;
	D1[7] = 3423500000.00;
	D1[8] = 2567140000.00;



	//set target values for Harris belly tuning: K = 4250, f5/f2 = 2.
	double target[2] = { 4250, 2 };
	
	//create thin solid
	BSplineSolid basis_solids[3];
	basis_solids[0] = solid;
	basis_solids[1] = solid;
	basis_solids[2] = solid;
	int index;
	double centre;
	double thin_factor = 0.5;
	for (int i = 0; i < solid.nx(); i++){
		for (int j = 0; j < solid.ny(); j++){
			centre = 0;
			for (int k = 0; k < solid.nz(); k++){
				index = i*solid.coffsetx() + j*solid.coffsety() + k * solid.coffsetz();
				centre += solid.cz[index]/solid.nz();
			}
			for (int k = 0; k < solid.nz(); k++){
				index = i*solid.coffsetx() + j*solid.coffsety() + k * solid.coffsetz();
				basis_solids[1].cz[index] = (centre + solid.cz[index])* thin_factor;
			}
		}
	}
	thin_factor = 0.5;
	for (int i = 0; i < solid.nx(); i++){
		if (i>solid.nx()*0.25 && i < solid.nx()*0.75)continue;
		for (int j = 0; j < solid.ny(); j++){
			if (j>solid.ny()*0.25 && j < solid.ny()*0.75)continue;
			centre = 0;
			for (int k = 0; k < solid.nz(); k++){
				index = i*solid.coffsetx() + j*solid.coffsety() + k * solid.coffsetz();
				centre += solid.cz[index] / solid.nz();
			}
			for (int k = 0; k < solid.nz(); k++){
				index = i*solid.coffsetx() + j*solid.coffsety() + k * solid.coffsetz();
				basis_solids[2].cz[index] = centre + (solid.cz[index] - centre)* thin_factor;
			}
		}
	}
	optimizer.tune_plate_harris(basis_solids, 3, density, D1, target);

}


bool TestOlsoAlgorithm(bool verbose){
	//create simple Bspline curve
	using namespace bspline;

	double tol = 0.0000001;

	double knot[6] = { 0.0, 0.0, 0.0, 1.0, 1.0, 1.0 };
	double cx[3] = { 0.0, 0.0, 1.0 };
	double cy[3] = { 0.0, 1.0, 1.0 };

	vec2<double> c[3];
	c[0] = vec2<double>(0.0, 0.0);
	c[1] = vec2<double>(0.0, 1.0);
	c[2] = vec2<double>(1.0, 1.0);
	int p = 2;
	int n = 3;

	int new_n = 10;
	double* new_knot = new double[new_n + p + 1];
	double* new_cx = new double[new_n];
	double* new_cy = new double[new_n];
	vec2<double>* new_c = new vec2<double>[new_n];
	for (int i = 0; i < new_n + p + 1; i++){
		double val = (double)(i - p) / (new_n - n + 1);
		if (val < 0){ new_knot[i] = 0; }
		else if (val > 1){ new_knot[i] = 1; }
		else{ new_knot[i] = val; }
	}

	bspline::olso_insertion(n - 1, cx, p + 1, knot, new_knot, new_n + p, new_cx);
	bspline::olso_insertion(n - 1, cy, p + 1, knot, new_knot, new_n + p, new_cy);
	bspline::olso_insertion(n - 1, c, p + 1, knot, new_knot, new_n + p, new_c);

	if (verbose){
		for (int i = 0; i < n; i++){
			cout << "old(cx, cy) = " << cx[i] << "," << cy[i] << endl;
		}
		for (int i = 0; i < new_n; i++){
			cout << "new(cx, cy) = " << new_cx[i] << "," << new_cy[i] << endl;
		}
		for (int i = 0; i < 20; i++){
			double x = (double)i / 19;
			double curve = 0;
			double new_curve = 0;
			for (int j = 0; j < n; j++){
				curve += cx[j] * bspline::deboor_value(j, p, knot, n + p + 1, x);
			}
			for (int j = 0; j < new_n; j++){
				new_curve += new_cx[j] * bspline::deboor_value(j, p, new_knot, new_n + p + 1, x);
			}
			cout << x << ": ";
			cout << bspline::evaluate_curve(p, *knot, n + p + 1, *cx, x) << ", ";
			cout << bspline::evaluate_curve(p, *new_knot, new_n + p + 1, *new_cx, x) << ", ";
			cout << curve << ", ";
			cout << new_curve << endl;
		}

	}

	for (int i = 0; i < 20; i++){
		double x = (double)i / 19;
		double px = bspline::evaluate_curve(p, *knot, n + p + 1, *cx, x);
		double py = bspline::evaluate_curve(p, *knot, n + p + 1, *cy, x);
		vec2<double> pv = bspline::evaluate_curve(p, *knot, n + p + 1, *c, x);

		double qx = bspline::evaluate_curve(p, *new_knot, new_n + p + 1, *new_cx, x);
		double qy = bspline::evaluate_curve(p, *new_knot, new_n + p + 1, *new_cy, x);
		vec2<double> qv = bspline::evaluate_curve(p, *new_knot, new_n + p + 1, *new_c, x);

		if (AreEqual(px, qx, tol)==false){ return false; }
		if (AreEqual(py, qy, tol) == false){ return false; }
		if (AreEqual(px, qv.x, tol) == false){ return false; }
		if (AreEqual(py, qv.y, tol) == false){ return false; }

	}


	//delete[] new_knot;
	delete[] new_cx;
	delete[] new_cy;
	delete[] new_c;



	//p = 1 spline
	double knot2[5] = { 0.0, 0.0, 0.5, 1.0, 1.0 };

	vec2<double> c2[3];
	c2[0] = vec2<double>(0.0, 0.0);
	c2[1] = vec2<double>(0.0, 1.0);
	c2[2] = vec2<double>(1.0, 1.0);
	p = 1;
	n = 3;
	 
	new_n = 7;
	double new_knot2[9] = { 0.0, 0.0, 0.25,0.25, 0.5, 0.75, 0.75, 1.0, 1.0 };
	vec2<double>* new_c2 = new vec2<double>[new_n];

	bspline::olso_insertion(n - 1, c2, p + 1, knot2, new_knot2, new_n + p, new_c2);

	if (verbose){
		for (int i = 0; i < n; i++){
			cout << "old(cx, cy) = " << c2[i].x << "," << c2[i].y << endl;
		}
		for (int i = 0; i < new_n; i++){
			cout << "new(cx, cy) = " << new_c2[i].x << "," << new_c2[i].y << endl;
		}
		for (int i = 0; i < 20; i++){
			double x = (double)i / 19;

			vec2<double> c_old = bspline::evaluate_curve(p, *knot2, n + p + 1, *c2, x);
			vec2<double> c_new = bspline::evaluate_curve(p, *new_knot2, new_n + p + 1, *new_c2, x);

			cout << x << ": " << c_old.x << ", " << c_old.y << ", " << c_new.x << ", " << c_new.y << endl;;
		}

	}

	for (int i = 0; i < 20; i++){
		double x = (double)i / 19;
		vec2<double> pv = bspline::evaluate_curve(p, *knot2, n + p + 1, *c2, x);
		vec2<double> qv = bspline::evaluate_curve(p, *new_knot2, new_n + p + 1, *new_c2, x);

		if (AreEqual(pv.x, qv.x, tol) == false){ return false; }
		if (AreEqual(pv.y, qv.y, tol) == false){ return false; }
	}
	delete[] new_c2;



	BSplineSolid input(1.0,1.0,1.0,2, 2, 2, 10, 5, 5);
	BSplineSolid output1;
	BSplineSolid output2;

	BSplineSolid::olso_split_spans_x(&input, &output1, 10);
	BSplineSolid::olso_split_spans_y(&output1, &output2, 10);
	
	if (verbose){
		cout << "input: degrees of freedom = ";
		cout << input.npoints() * 3 << endl;

		cout << "output1: degrees of freedom = ";
		cout << output1.npoints() * 3 << endl;

		cout << "output2: degrees of freedom = ";
		cout << output2.npoints() * 3 << endl;
	}

	bspline::curve<double> curve(2,0.5);
	curve.append(1.0);
	curve.append(2.5);
	double x = curve.evaluate(0.5);

	bspline::vec2<float> v1(0.5f, 0.5f);
	bspline::vec2<float> v2(0.1f, 1.0f);
	bspline::vec2<float> v3(1.0f, 2.5f);
	bspline::curve<bspline::vec2<float>> curve2(2, v1);
	curve2.append(v2);
	curve2.append(v3);

	bspline::vec2<float> v4 = curve2.evaluate(0.5);

	//bspline::vec2<float> res1 = curve2.evaluate(0.25);
	//bspline::vec2<float> res2 = curve2.evaluate(0.5);
	//bspline::vec2<float> res3 = curve2.evaluate(0.75);

	return true;
}



bool TestSurface(bool verbose)
{

	SURFACE3F s1(
		2,
		1,
		std::vector<double>{ 0, 0, 0, 1, 1, 1 },
		std::vector<double>{ 0, 0, 1, 1 },
		LATTICE3F(std::vector<VEC3F>{ VEC3F(0, 0, 0), VEC3F(0.3, 0, 0), VEC3F(1.0, 0, 0), VEC3F(0, 1.0, 0), VEC3F(0.3, 1.0, 0), VEC3F(1.0, 1.0, 0) }, 3, 2)
	);
	
	SURFACE3F s2(s1);
	s2.p = 4;
	s2.knotx[0] = -1;
	s2.points.data[0].x = -2;

	VEC3F p1 = s1.evaluate(0, 0);
	VEC3F p2 = s1.evaluate(1, 1);
	VEC3F p3 = s1.evaluate(0.5, 0.5);


	if (true){
		return true;
	}
	else{
		return false;
	}

}


bool TestCurve_InsertKnots(bool verbose){

	using namespace bspline;

	int p = 2;
	int n = 5;
	float r = 1;
	float theta = 0;
	curve<vec2<float>>* c = new curve<vec2<float>>(p, vec2<float>(r*cosf(theta), r*sinf(theta)));
	for (int i = 0; i < n - 1; i++){
		theta += 0.1;
		c->append(vec2<float>(r*cosf(theta), r*sinf(theta)));
	}

	double newKnots[] = { 0.1, 0.9 };
	curve<vec2<float>>* d;
	d = c->insert_knots(newKnots, 2);
	
	curve<vec2<float>>* e;
	e = c->trim_curve(newKnots[0], newKnots[1]);


	int neval = 15;
	double t;

	if (verbose){
		for (int i = 0; i <= neval; i++){
			t = (i * c->maxParam() + neval * c->minParam()) / neval;
			cout << t << ": ";
			cout << "(" << c->evaluate(t).x << "," << c->evaluate(t).y << ") ";
			cout << "(" << d->evaluate(t).x << "," << d->evaluate(t).y << ") ";
			cout << endl;
		}


		for (int i = 0; i < c->lKnot(); i++){
			cout << "c: " << c->getKnot(i) << endl;
		}
		cout << endl;
		for (int i = 0; i < d->lKnot(); i++){
			cout << "d: " << d->getKnot(i) << endl;
		}

		for (int i = 0; i < c->nPoints(); i++){
			cout << "c[" << i << "] = (" << c->get(i).x << "," << c->get(i).y << ") " << endl;
		}
		for (int i = 0; i < d->nPoints(); i++){
			cout << "d[" << i << "] = (" << d->get(i).x << "," << d->get(i).y << ") " << endl;
		}

		cout << "Trimmed curve" << endl;
		cout << "Original:  New:" << endl;
		for (int i = 0; i <= neval; i++){
			t = (i * e->maxParam() + neval * e->minParam()) / neval;
			cout << t << ": ";
			cout << "(" << c->evaluate(t).x << "," << c->evaluate(t).y << ") ";
			cout << "(" << e->evaluate(t).x << "," << e->evaluate(t).y << ") ";
			cout << endl;
		}
	}

	double tol = 0.000001;
	if (c->maxParam() != d->maxParam()){ return false; }
	if (c->minParam() != d->minParam()){ return false; }
	if (e->minParam() != newKnots[0]){ return false; }
	if (e->maxParam() != newKnots[1]){ return false; }
	if (c->nPoints() != d->nPoints() - 2){ return false; }
	if (c->lKnot() != d->lKnot() - 2){ return false; }
	for (int i = 0; i <= neval; i++){
		t = (i * c->maxParam() + neval * c->minParam()) / neval;
		vec2<float> p = c->evaluate(t);
		vec2<float> q = d->evaluate(t);
		if (AreEqual(p.x, q.x, tol)==false){ return false; }
		if (AreEqual(p.y, q.y, tol) == false){ return false; }
	}
	for (int i = 0; i <= neval; i++){
		t = (i * e->maxParam() + neval * e->minParam()) / neval;
		vec2<float> p = c->evaluate(t);
		vec2<float> q = e->evaluate(t);
		if (AreEqual(p.x, q.x, tol) == false){ return false; }
		if (AreEqual(p.y, q.y, tol) == false){ return false; }
	}
	return true;
}


bool TestCurve_LineIntersection(bool verbose){

	using namespace bspline;

	vec2<float> p1(0,0);
	vec2<float> p2(10, 10);
	vec2<float> q1(5, 2);
	vec2<float> q2(5, 6);

	vec2<double> params;
	bool intersect;

	intersect = bspline::IntersectionParam(p1,p2,q1,q2,&params,verbose);
	if (intersect == false || params.x != 0.5 || params.y != 0.25){ return false; }

	intersect = bspline::IntersectionParam(q1, q2, p1, p2, &params, verbose);
	if (intersect == false || params.x != 0.25 || params.y != 0.5){ return false; }

	q1.y = 5;
	intersect = bspline::IntersectionParam(p1, p2, q1, q2, &params, verbose);
	if (intersect == false || params.x != 0.5 || params.y != 1.0){ return false; }

	q1.y = 0;
	q2.y = 5;
	intersect = bspline::IntersectionParam(p1, p2, q1, q2, &params, verbose);
	if (intersect == false || params.x != 0.5 || params.y != 0.0){ return false; }

	q2.y = 4;
	intersect = bspline::IntersectionParam(p1, p2, q1, q2, &params, verbose);
	if (intersect == true ){ return false; }

	q1.x = 5;
	q1.y = 0;
	q2.x = 10;
	q2.y = 5;
	intersect = bspline::IntersectionParam(p1, p2, q1, q2, &params, verbose);
	if (intersect == true){ return false; }

	q1.x = 5;
	q1.y = 5;
	q2.x = 7;
	q2.y = 7;
	intersect = bspline::IntersectionParam(p1, p2, q1, q2, &params, verbose);
	if (intersect == true){ return false; }

	return true;
}


bool TestCurve_CurveIntersection(bool verbose){

	using namespace bspline;

	int p = 2;
	int n = 5;
	float r = 1;
	float theta = 0;
	curve<vec2<float>>* c1 = new curve<vec2<float>>(p, vec2<float>(r*cosf(theta), r*sinf(theta)));
	curve<vec2<float>>* c2 = new curve<vec2<float>>(p, vec2<float>(r - r*cosf(-theta), r*sinf(theta)));
	for (int i = 0; i < n - 1; i++){
		theta += 0.6;
		c1->append(vec2<float>(r*cosf(theta), r*sinf(theta)));
		c2->append(vec2<float>(r - r*cosf(-theta), r*sinf(theta)));
	}

	if (verbose){
		for (int i = 0; i < c1->nPoints(); i++){
			cout << "i = " << i;
			cout << " c1 = (" << c1->get(i).x << "," << c1->get(i).y;
			cout << " c2 = (" << c2->get(i).x << "," << c2->get(i).y;
			cout << endl;
		}
	}

	vec2<double> param;
	float error = 0.001;
	bool intersectFound = IntersectParam((*c1), (*c2), error, &param, verbose);

	vec2<float> i1 = c1->evaluate(param.x);
	vec2<float> i2 = c2->evaluate(param.y);
	float d = bspline::distance(i1, i2);
	if (intersectFound == false){ return false; }
	if (d > error){ return false; }
	return true;

}


bool TestTriangle(bool verbose){

	using namespace bspline;

	VEC2F p1(1.0, 0.0);
	VEC2F p2(9.0, 0.0);
	VEC2F p3(1.0, 7.0);

	VEC2F c1(0.0, 0.0);
	float r1 = 0.9;

	bool overlap1 = interior_triangle(c1, r1, p1, p2, p3);

	VEC2F c2(0.0, 0.0);
	float r2 = 1.0;

	bool overlap2 = interior_triangle(c2, r2, p1, p2, p3);

	VEC2F c3(0.0, 0.0);
	float r3 = 1.1;

	bool overlap3 = interior_triangle(c3, r3, p1, p2, p3);

	VEC2F c4(3.0, 1.5);
	float r4 = 1.0;

	bool overlap4 = interior_triangle(c4, r4, p1, p2, p3);

	if (overlap1 == true){ return false; }
	if (overlap2 == false){ return false; }
	if (overlap3 == false){ return false; }
	if (overlap4 == false){ return false; }
	return true;

}




bool TestToolPath3(bool verbose){

	char filepath[] = "C:\\Users\\Lizzie\\Documents\\GitHub\\Finite-Element-Engine\\FiniteElementAnalysis\\Data\\Perlman Strad Violin Model.xml";

	app_model* model = IO::LoadModelXML(filepath, verbose);
	if (model == nullptr){ return false; }
	if (model->violin == nullptr){ return false; }

	model->violin->ribs->rib_mould_locator_holes.push_back(circle2f(-140, -65, 10));
	model->violin->ribs->rib_mould_locator_holes.push_back(circle2f(-140, 65, 10));
	model->violin->ribs->rib_mould_locator_holes.push_back(circle2f(110, 45, 10));
	model->violin->ribs->rib_mould_locator_holes.push_back(circle2f(110, -45, 10));

	for (std::pair<std::string, toolpath_base*> path : model->paths){
		path.second->calculate();
		path.second->save_gcode();
	}

	return true;
}

