#include "arpack_shiftinvert.h"
#include <stdlib.h>   
#include <malloc.h>
#include <memory.h>
#include <iostream>
#include <fstream>
#include "slu_ddefs.h"
#include <math.h>

arpack_shiftinvert::arpack_shiftinvert(SuperMatrix &inputM, SuperMatrix &inputA)
{
	
	A = inputA;
	M = inputM;
	
	//main options:
	SIGMA = 0;
	NEV = 20;
	RVEC = true;
	HOWMNY = 'A';
	WHICH[0] = 'L';
	WHICH[1] = 'M';
	int ishfts = 1;
	int	maxitr = 300;
	int mode = 3;

	//note fortran indeces start at 1, C++ is zero-based
	IPARAM[0] = ishfts;
	IPARAM[2] = maxitr;
	IPARAM[6] = mode;


	N = A.nrow;
	NCV = 4 * NEV;
	LWORKL = NCV * (NCV + 8);
	LDV = N;
	LDZ = N;

	int sizeofbool = sizeof(bool);
	int sizeofint = sizeof(int);

	SELECT = (int*)malloc(NCV * sizeof(int));
	RESID = (double*)malloc(N * sizeof(double));
	for (int i = 0; i < N; i++){
		RESID[i] = 0;
	}
	WORKD = (double*)malloc(3 * N * sizeof(double));
	WORKL = (double*)malloc(LWORKL * sizeof(double));
	V = (double*)malloc(N * NCV * sizeof(double));
	D = (double*)malloc(NEV * sizeof(double));
	
	if (RVEC){
		Z = (double*)malloc(N * NEV * sizeof(double));
	}
	

	//SuperLU objects
	m = A.nrow;
	n = A.ncol;
	
	/* Defaults */
	lwork = 0;
	nrhs = 1;
	equil = YES;
	u = 1.0;
	trans = NOTRANS;
	
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
	options.IterRefine = NOREFINE;//SLU_DOUBLE;

	if (lwork > 0) {
		work = SUPERLU_MALLOC(lwork);
		if (!work) {
			ABORT("DLINSOLX: cannot allocate work[]");
		}
	}

	if (!(rhsb = doubleMalloc(m * nrhs))) ABORT("Malloc fails for rhsb[].");
	if (!(rhsx = doubleMalloc(m * nrhs))) ABORT("Malloc fails for rhsx[].");
	dCreate_Dense_Matrix(&B, m, nrhs, rhsb, m, SLU_DN, SLU_D, SLU_GE);
	dCreate_Dense_Matrix(&X, m, nrhs, rhsx, m, SLU_DN, SLU_D, SLU_GE);
	xact = doubleMalloc(n * nrhs);
	ldx = n;
	dGenXtrue(n, nrhs, xact, ldx);
	dFillRHS(trans, nrhs, xact, ldx, &A, &B);

	if (!(etree = intMalloc(n))) ABORT("Malloc fails for etree[].");
	if (!(perm_r = intMalloc(m))) ABORT("Malloc fails for perm_r[].");
	if (!(perm_c = intMalloc(n))) ABORT("Malloc fails for perm_c[].");
	if (!(R = (double *)SUPERLU_MALLOC(A.nrow * sizeof(double))))
		ABORT("SUPERLU_MALLOC fails for R[].");
	if (!(C = (double *)SUPERLU_MALLOC(A.ncol * sizeof(double))))
		ABORT("SUPERLU_MALLOC fails for C[].");
	if (!(ferr = (double *)SUPERLU_MALLOC(nrhs * sizeof(double))))
		ABORT("SUPERLU_MALLOC fails for ferr[].");
	if (!(berr = (double *)SUPERLU_MALLOC(nrhs * sizeof(double))))
		ABORT("SUPERLU_MALLOC fails for berr[].");

	/* Initialize the statistics variables. */
	StatInit(&stat);


}


arpack_shiftinvert::~arpack_shiftinvert()
{
	//ARPACK arrays
	free(SELECT);
	free(RESID);
	free(WORKD);
	free(WORKL);
	free(V);
	free(D);
	if (RVEC){
		free(Z);
	}


	//superLU arrays
	SUPERLU_FREE(rhsb);
	SUPERLU_FREE(rhsx);
	SUPERLU_FREE(xact);
	SUPERLU_FREE(etree);
	SUPERLU_FREE(perm_r);
	SUPERLU_FREE(perm_c);
	SUPERLU_FREE(R);
	SUPERLU_FREE(C);
	SUPERLU_FREE(ferr);
	SUPERLU_FREE(berr);
	Destroy_CompCol_Matrix(&A);
	Destroy_SuperMatrix_Store(&B);
	Destroy_SuperMatrix_Store(&X);
	if (lwork == 0) {
		Destroy_SuperNode_Matrix(&L);
		Destroy_CompCol_Matrix(&U);
	}
	else if (lwork > 0) {
		SUPERLU_FREE(work);
	}

}


void arpack_shiftinvert::setM(SuperMatrix& inputM)
{
	if (inputM.Dtype != M.Dtype)ABORT("setM fails: incompatible Dtype."); 
	if (inputM.Mtype != M.Mtype)ABORT("setM fails: incompatible Mtype.");
	if (inputM.Stype != M.Stype)ABORT("setM fails: incompatible Stype.");
	if (inputM.ncol != M.ncol)ABORT("setM fails: incompatible ncol.");
	if (inputM.nrow != M.nrow)ABORT("setM fails: incompatible nrow.");
	
	M = inputM;
}
void arpack_shiftinvert::setA(SuperMatrix& inputA)
{
	if (inputA.Dtype != A.Dtype)ABORT("setA fails: incompatible Dtype.");
	if (inputA.Mtype != A.Mtype)ABORT("setA fails: incompatible Mtype.");
	if (inputA.Stype != A.Stype)ABORT("setA fails: incompatible Stype.");
	if (inputA.ncol != A.ncol)ABORT("setA fails: incompatible ncol.");
	if (inputA.nrow != A.nrow)ABORT("setA fails: incompatible nrow.");

	A = inputA;
}



void arpack_shiftinvert::solve()
{
	int count = 0;
	int IDO = 0;
	INFO = 0;


	bool continueLoop = true;
	while(continueLoop){

		//call arpack function for shift-inverted lanczos algorithm
		//print(IDO);
		
		ARPACK_dsaupd(&IDO, &BMAT, &N, WHICH, &NEV, &TOL, RESID, &NCV, V, &LDV, IPARAM, IPNTR, WORKD, WORKL, &LWORKL, &INFO);

		switch (IDO)
		{
			
		case -1:
			/*
			Initialisation: Perform  y <-- - OP * x = inv[A - SIGMA * M] * M * x
			to force the starting vector into the range of OP. The matrix vector multiplication routine
			takes workd(ipntr(0)) as the input vector. 
			The final result is returned to            
			workd(ipntr(1))
			*/
			multiply(M, &WORKD[IPNTR[0] - 1], &WORKD[IPNTR[1] - 1]);
			Bstore = (DNformat *)B.Store;
			Bstore->nzval = &WORKD[IPNTR[1] - 1];
			
			/* ------------------------------------------------------------
			NOW WE SOLVE THE LINEAR SYSTEM USING THE FACTORED FORM OF A.
			------------------------------------------------------------*/
			if (!availableLUfact) factorizeLU();
			
			options.Fact = FACTORED; /* Indicate the factored form of A is supplied. */
			B.ncol = nrhs;  /* Set the number of right-hand side */

			/* Initialize the statistics variables. */
			StatInit(&stat);

			dgssvx(&options, &A, perm_c, perm_r, etree, equed, R, C,
				&L, &U, work, lwork, &B, &X, &rpg, &rcond, ferr, berr,
				&Glu, &mem_usage, &stat, &info);

			memcpy(&WORKD[IPNTR[1] - 1], (double*)((DNformat*)X.Store)->nzval, N * sizeof(double));
			break;
		case 1:
			/*
			Perform y <-- OP*x = inv[A - sigma * M] * M * x 
			M*x has been saved in workd(ipntr(2)). The user only needs the linear system   
			solver here that takes workd(ipntr(2) as input, and returns the result to     
			workd(ipntr(1)).
			*/
			Bstore->nzval = &WORKD[IPNTR[2] - 1];
			StatInit(&stat);
			dgssvx(&options, &A, perm_c, perm_r, etree, equed, R, C,
				&L, &U, work, lwork, &B, &X, &rpg, &rcond, ferr, berr,
				&Glu, &mem_usage, &stat, &info);
			
			memcpy(&WORKD[IPNTR[1] - 1], (double*)((DNformat*)X.Store)->nzval, N * sizeof(double));
			if (info == 0) {

			}
			break;
		case 2:
			/*
			Perform  y <--- M*x
			Need the matrix vector multiplication routine here that takes workd(ipntr(0))
			as the input and returns the result to workd(ipntr(1)).
			*/
			multiply(M, &WORKD[IPNTR[0] - 1], &WORKD[IPNTR[1] - 1]);
			break;
		default:
			continueLoop = false;
		}
	} 

	ARPACK_dseupd(&RVEC, &HOWMNY, SELECT, D, Z, &LDZ, &SIGMA, &BMAT, &N, WHICH, &NEV, &TOL, RESID, &NCV, V, &LDV, IPARAM, IPNTR, WORKD, WORKL, &LWORKL, &INFO);

	std::cout << "Info = " << INFO << "\n";
	std::cout << "Number of converged eigenvalues = " << IPARAM[4] << "\n";

	std::cout << "Eigenvalues:\n";
	for (int i = 0; i < IPARAM[4]; i++){
		std::cout << sqrt(D[i]) / 3.141592653589793238 / 2 << "\n";
	}

}


double arpack_shiftinvert::get_frequency(int i){
	if (i < 1)return 0;
	if (i>=IPARAM[4] - 5)return 0;
	return sqrt(D[i+5]) / 3.141592653589793238 / 2;
}


void arpack_shiftinvert::factorizeLU(){
	
	/* ONLY PERFORM THE LU DECOMPOSITION */
	B.ncol = 0;  /* Indicate not to solve the system */
	options.Fact = DOFACT;
	dgssvx(&options, &A, perm_c, perm_r, etree, equed, R, C,
		&L, &U, work, lwork, &B, &X, &rpg, &rcond, ferr, berr,
		&Glu, &mem_usage, &stat, &info);

	printf("LU factorization: dgssvx() returns info %d\n", info);

	if (info == 0 || info == n + 1) {

		if (options.PivotGrowth) printf("Recip. pivot growth = %e\n", rpg);
		if (options.ConditionNumber)
			printf("Recip. condition number = %e\n", rcond);
		Lstore = (SCformat *)L.Store;
		Ustore = (NCformat *)U.Store;
		printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
		printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
		printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz - n);
		printf("FILL ratio = %.1f\n", (float)(Lstore->nnz + Ustore->nnz - n) / nnz);

		printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
			mem_usage.for_lu / 1e6, mem_usage.total_needed / 1e6);
		fflush(stdout);

		availableLUfact = true;
	}
	else if (info > 0 && lwork == -1) {
		printf("** Estimated memory: %d bytes\n", info - n);
	}

	if (options.PrintStat) StatPrint(&stat);
	StatFree(&stat);

}


void arpack_shiftinvert::multiply(SuperMatrix &M, const double* x, double* y){
	NCformat *Mstore = (NCformat*)M.Store;
	double *nzval = (double*)Mstore->nzval;
	for (int i = 0; i < M.nrow; i++){
		y[i] = 0;
		for (int j = Mstore->colptr[i]; j < Mstore->colptr[i + 1]; j++){
			y[i] += nzval[j] * x[Mstore->rowind[j]];
		}
	}
}

void arpack_shiftinvert::shift(const double sigma){
	SIGMA = sigma;
	NCformat *Mstore = (NCformat*)M.Store;
	double *Mnzval = (double*)Mstore->nzval;

	NCformat *Astore = (NCformat*)A.Store;
	double *Anzval = (double*)Astore->nzval;

	for (int i = 0; i < Astore->nnz; i++){
		Anzval[i] -= sigma * Mnzval[i];
	}
}


void arpack_shiftinvert::print(int IDO){

	std::ofstream outfile;
	outfile.open("output.txt", std::ios_base::app);

	outfile << "IDO = " << IDO << "\n";
	outfile << "SIGMA = " << SIGMA << "\n";
	outfile << "NEV = " << NEV << "\n";
	outfile << "RVEC = " << RVEC << "\n";
	outfile << "HOWMNY = " << HOWMNY << "\n";
	outfile << "WHICH = " << WHICH[0] << WHICH[1] << "\n";

	outfile << "IPARAM = ";
	for (int i = 0; i < 11; i++){
		outfile << IPARAM[i] << ",";
	}
	outfile << "\n";

	outfile << "N = " << N << "\n";
	outfile << "NCV = " << NCV << "\n";
	outfile << "LWORKL = " << LWORKL << "\n";
	outfile << "LDZ = " << LDZ << "\n";

	outfile << "SELECT = ";
	for (int i = 0; i < NCV; i++){
		outfile << SELECT[i] << ",";
	}
	outfile << "\n";

	outfile << "RESID = ";
	for (int i = 0; i < N; i++){
		outfile << RESID[i] << ",";
	}
	outfile << "\n";

	outfile << "WORKD = ";
	for (int i = 0; i < 3*N; i++){
		outfile << WORKD[i] << ",";
	}
	outfile << "\n";

	outfile << "WORKL = ";
	for (int i = 0; i < LWORKL; i++){
		outfile << WORKL[i] << ",";
	}
	outfile << "\n";

	outfile << "V = ";
	for (int i = 0; i < N*NCV; i++){
		outfile << V[i] << ",";
	}
	outfile << "\n";

	outfile << "D = ";
	for (int i = 0; i < NEV; i++){
		outfile << D[i] << ",";
	}
	outfile << "\n";

	if (RVEC){
		outfile << "Z = ";
		for (int i = 0; i < N * NEV; i++){
			outfile << Z[i] << ",";
		}
		outfile << "\n";
	}

	outfile << "\n\n\n";



}