#pragma once

#include "slu_ddefs.h"


/*
Fortran functions take all parameters as pointers, so the C++ function must do the same
extern "C" {....} is required so that C++ compiler doesn't change the function names
the gfortran compiler appended an underscore to the function name (name-mangling), so
C function must have underscore appended, i.e. function dsaupd in fortran code referred to as dsaupd_ in C.

For .exe to run, it needs x64 versions of the mingw .dlls in the same folder
*/

#ifndef ARPACK_HEADER_INCLUDED
#define ARPACK_HEADER_INCLUDED

#ifndef ARPACK_GLOBAL
#if defined(ARPACK_GLOBAL_PATTERN_LC) || defined(ADD_)
#define ARPACK_GLOBAL(lcname,UCNAME)  lcname##_
#elif defined(ARPACK_GLOBAL_PATTERN_UC) || defined(UPPER)
#define ARPACK_GLOBAL(lcname,UCNAME)  UCNAME
#elif defined(ARPACK_GLOBAL_PATTERN_MC) || defined(NOCHANGE)
#define ARPACK_GLOBAL(lcname,UCNAME)  lcname
#else
#define ARPACK_GLOBAL(lcname,UCNAME)  lcname##_
#endif
#endif

#endif

extern "C" {

#define ARPACK_dsaupd ARPACK_GLOBAL(dsaupd,DSAUPD)
#define ARPACK_dseupd ARPACK_GLOBAL(dseupd,DSEUPD)

	//void __cdecl ARPACK_dsaupd(int* IDO, char* BMAT, int* N, char* WHICH, int* NEV, double* TOL, double* RESID, int* NCV, double* V, int* LDV, int* IPARAM, int* IPNTR, double* WORKD, double* WORKL, int* LWORKL, int* INFO);
	//void __cdecl ARPACK_dseupd(bool *RVEC, char *HOWMNY, int *SELECT, double *D, double *Z, int *LDZ, double *SIGMA, char *BMAT, int *N, char *WHICH, int *NEV, double *TOL, double* RESID, int *NCV, double *V, int* LDV, int *IPARAM, int *IPNTR, double *WORKD, double *WORKL, int *LWORKL, int *INFO);

	void __cdecl ARPACK_dsaupd(int* IDO, char* BMAT, int* N, char* WHICH, int* NEV, double* TOL, double* RESID, int* NCV, double* V, int* LDV, int* IPARAM, int* IPNTR, double* WORKD, double* WORKL, int* LWORKL, int* INFO);
	void __cdecl ARPACK_dseupd(bool *RVEC, char *HOWMNY, int *SELECT, double *D, double *Z, int *LDZ, double *SIGMA, char *BMAT, int *N, char *WHICH, int *NEV, double *TOL, double* RESID, int *NCV, double *V, int* LDV, int *IPARAM, int *IPNTR, double *WORKD, double *WORKL, int *LWORKL, int *INFO);

}
//
//#ifdef FEADLL_EXPORTS
//#define FEADLL_API __declspec(dllexport) 
//#else
//#define FEADLL_API __declspec(dllimport) 
//#endif

class arpack_shiftinvert
{
public:
	arpack_shiftinvert(SuperMatrix&, SuperMatrix&);
	~arpack_shiftinvert();

	void arpack_shiftinvert::solve();
	void arpack_shiftinvert::shift(double);
	void arpack_shiftinvert::factorizeLU();
	void arpack_shiftinvert::setM(SuperMatrix&);
	void arpack_shiftinvert::setA(SuperMatrix&);

	/*
	Returns frequency of ith mode, indexed 1, 2, 3 ... N
	*/
	double arpack_shiftinvert::get_frequency(int i);

	SuperMatrix M;
	SuperMatrix A;

private:
	
	void arpack_shiftinvert::multiply(SuperMatrix&, const double*, double*);
	void arpack_shiftinvert::print(int IDO);

	bool availableLUfact = false;

	//SuperLU objects
	char           equed[1];
	yes_no_t       equil;
	trans_t        trans;
	SuperMatrix    L, U;
	SuperMatrix    B, X;
	NCformat       *Astore;
	NCformat       *Ustore;
	DNformat       *Bstore;
	DNformat       *Xstore;
	SCformat       *Lstore;
	GlobalLU_t	   Glu; /* facilitate multiple factorizations with
						SamePattern_SameRowPerm                  */
	double         *a;
	int            *asub, *xa;
	int            *perm_c; /* column permutation vector */
	int            *perm_r; /* row permutations from partial pivoting */
	int            *etree;
	void           *work;
	int            info, lwork, nrhs, ldx;
	int            i, m, n, nnz;
	double         *rhsb, *rhsx, *xact;
	double         *R, *C;
	double         *ferr, *berr;
	double         u, rpg, rcond;
	mem_usage_t    mem_usage;
	superlu_options_t options;
	SuperLUStat_t stat;

	/*
	BMAT specifies the type of the matrix B that defines the
	semi - inner product for the operator OP.
	B = 'I'->standard eigenvalue problem A*x = lambda*x
	B = 'G'->generalized eigenvalue problem A*x = lambda*B*x
	*/
	char BMAT = 'G';

	//dimension of eigenproblem
	int N;	


	/*
	Specify which of the Ritz values of OP to compute.

    'LA' - compute the NEV largest (algebraic) eigenvalues.
    'SA' - compute the NEV smallest (algebraic) eigenvalues.
    'LM' - compute the NEV largest (in magnitude) eigenvalues.
	'SM' - compute the NEV smallest (in magnitude) eigenvalues. 
	'BE' - compute NEV eigenvalues, half from each end of the
			spectrum.  When NEV is odd, compute one more from the
			high end than from the low end.
	*/
	char WHICH[2];		

	//Number of eigenvalues of OP to be computed. 0 < NEV < N.
	int NEV;


	/*
	Stopping criterion: the relative accuracy of the Ritz value 
	is considered acceptable if BOUNDS(I) .LE. TOL*ABS(RITZ(I)).
	If TOL <= 0. is passed a default is set:
	DEFAULT = DLAMCH ('EPS')  (machine precision as computed
	by the LAPACK auxiliary subroutine DLAMCH ).
	*/
	double TOL = 0;


	/*
	Double precision  array of length N.  (INPUT/OUTPUT)
	On INPUT: 
	If INFO .EQ. 0, a random initial residual vector is used.
	If INFO .NE. 0, RESID contains the initial residual vector,possibly from a previous run.
	On OUTPUT:
	RESID contains the final residual vector. 
	*/
	double *RESID;

	/*
	Number of columns of the matrix V (less than or equal to N).
c          This will indicate how many Lanczos vectors are generated 
c          at each iteration.  After the startup phase in which NEV 
c          Lanczos vectors are generated, the algorithm generates 
c          NCV-NEV Lanczos vectors at each subsequent update iteration.
c          Most of the cost in generating each Lanczos vector is in the 
c          matrix-vector product OP*x
	*/
	int NCV;

	/*
	Double precision  N by NCV array.  (OUTPUT)
	The NCV columns of V contain the Lanczos basis vectors.
	*/
	double *V;


	/*
	Leading dimension of V exactly as declared in the calling program.
	*/
	int LDV;


	/*
	Integer array of length 11.  (INPUT/OUTPUT)
c          IPARAM(1) = ISHIFT: method for selecting the implicit shifts.
c          The shifts selected at each iteration are used to restart
c          the Arnoldi iteration in an implicit fashion.
c          -------------------------------------------------------------
c          ISHIFT = 0: the shifts are provided by the user via
c                      reverse communication.  The NCV eigenvalues of
c                      the current tridiagonal matrix T are returned in
c                      the part of WORKL array corresponding to RITZ.
c                      See remark 6 below.
c          ISHIFT = 1: exact shifts with respect to the reduced 
c                      tridiagonal matrix T.  This is equivalent to 
c                      restarting the iteration with a starting vector 
c                      that is a linear combination of Ritz vectors 
c                      associated with the "wanted" Ritz values.
c          -------------------------------------------------------------
c
c          IPARAM(2) = LEVEC
c          No longer referenced. See remark 2 below.
c
c          IPARAM(3) = MXITER
c          On INPUT:  maximum number of Arnoldi update iterations allowed. 
c          On OUTPUT: actual number of Arnoldi update iterations taken. 
c
c          IPARAM(4) = NB: blocksize to be used in the recurrence.
c          The code currently works only for NB = 1.
c
c          IPARAM(5) = NCONV: number of "converged" Ritz values.
c          This represents the number of Ritz values that satisfy
c          the convergence criterion.
c
c          IPARAM(6) = IUPD
c          No longer referenced. Implicit restarting is ALWAYS used. 
c
c          IPARAM(7) = MODE
c          On INPUT determines what type of eigenproblem is being solved.
c          Must be 1,2,3,4,5; See under \Description of dsaupd  for the 
c          five modes available.
c
c          IPARAM(8) = NP
c          When ido = 3 and the user provides shifts through reverse
c          communication (IPARAM(1)=0), dsaupd  returns NP, the number
c          of shifts the user is to provide. 0 < NP <=NCV-NEV. See Remark
c          6 below.
c
c          IPARAM(9) = NUMOP, IPARAM(10) = NUMOPB, IPARAM(11) = NUMREO,
c          OUTPUT: NUMOP  = total number of OP*x operations,
c                  NUMOPB = total number of B*x operations if BMAT='G',
c                  NUMREO = total number of steps of re-orthogonalization. 
	*/
	int IPARAM[11];


	/*
	Integer array of length 11.  (OUTPUT)
c          Pointer to mark the starting locations in the WORKD and WORKL
c          arrays for matrices/vectors used by the Lanczos iteration.
c          -------------------------------------------------------------
c          IPNTR(1): pointer to the current operand vector X in WORKD.
c          IPNTR(2): pointer to the current result vector Y in WORKD.
c          IPNTR(3): pointer to the vector B * X in WORKD when used in 
c                    the shift-and-invert mode.
c          IPNTR(4): pointer to the next available location in WORKL
c                    that is untouched by the program.
c          IPNTR(5): pointer to the NCV by 2 tridiagonal matrix T in WORKL.
c          IPNTR(6): pointer to the NCV RITZ values array in WORKL.
c          IPNTR(7): pointer to the Ritz estimates in array WORKL associated
c                    with the Ritz values located in RITZ in WORKL.
c          IPNTR(11): pointer to the NP shifts in WORKL. See Remark 6 below.
c
c          Note: IPNTR(8:10) is only referenced by dseupd . See Remark 2.
c          IPNTR(8): pointer to the NCV RITZ values of the original system.
c          IPNTR(9): pointer to the NCV corresponding error bounds.
c          IPNTR(10): pointer to the NCV by NCV matrix of eigenvectors
c                     of the tridiagonal matrix T. Only referenced by
c                     dseupd  if RVEC = .TRUE. See Remarks.
	*/
	int IPNTR[11];



	/*
	Double precision  work array of length 3*N.  (REVERSE COMMUNICATION)
c          Distributed array to be used in the basic Arnoldi iteration
c          for reverse communication.  The user should not use WORKD 
c          as temporary workspace during the iteration. Upon termination
c          WORKD(1:N) contains B*RESID(1:N). If the Ritz vectors are desired
c          subroutine dseupd  uses this output.
c          See Data Distribution Note below.
	*/
	double *WORKD;


	/*
	Double precision  work array of length LWORKL.  (OUTPUT/WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.  See Data Distribution Note below.
	*/
	double *WORKL;


	/*
	LWORKL must be at least NCV**2 + 8*NCV 
	*/
	int LWORKL;


	/*
c	DSAUPD: INFO    Integer.  (INPUT/OUTPUT)
c          If INFO .EQ. 0, a randomly initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          Error flag on output.
c          =  0: Normal exit.
c          =  1: Maximum number of iterations taken.
c                All possible eigenvalues of OP has been found. IPARAM(5)  
c                returns the number of wanted converged Ritz values.
c          =  2: No longer an informational error. Deprecated starting
c                with release 2 of ARPACK.
c          =  3: No shifts could be applied during a cycle of the 
c                Implicitly restarted Arnoldi iteration. One possibility 
c                is to increase the size of NCV relative to NEV. 
c                See remark 4 below.
c          = -1: N must be positive.
c          = -2: NEV must be positive.
c          = -3: NCV must be greater than NEV and less than or equal to N.
c          = -4: The maximum number of Arnoldi update iterations allowed
c                must be greater than zero.
c          = -5: WHICH must be one of 'LM', 'SM', 'LA', 'SA' or 'BE'.
c          = -6: BMAT must be one of 'I' or 'G'.
c          = -7: Length of private work array WORKL is not sufficient.
c          = -8: Error return from trid. eigenvalue calculation;
c                Informatinal error from LAPACK routine dsteqr .
c          = -9: Starting vector is zero.
c          = -10: IPARAM(7) must be 1,2,3,4,5.
c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatable.
c          = -12: IPARAM(1) must be equal to 0 or 1.
c          = -13: NEV and WHICH = 'BE' are incompatable.
c          = -9999: Could not build an Arnoldi factorization.
c                   IPARAM(5) returns the size of the current Arnoldi
c                   factorization. The user is advised to check that
c                   enough workspace and array storage has been allocated.


c  DSEUPD: INFO    Integer.  (OUTPUT)
c          Error flag on output.
c          =  0: Normal exit.
c          = -1: N must be positive.
c          = -2: NEV must be positive.
c          = -3: NCV must be greater than NEV and less than or equal to N.
c          = -5: WHICH must be one of 'LM', 'SM', 'LA', 'SA' or 'BE'.
c          = -6: BMAT must be one of 'I' or 'G'.
c          = -7: Length of private work WORKL array is not sufficient.
c          = -8: Error return from trid. eigenvalue calculation;
c                Information error from LAPACK routine dsteqr .
c          = -9: Starting vector is zero.
c          = -10: IPARAM(7) must be 1,2,3,4,5.
c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
c          = -12: NEV and WHICH = 'BE' are incompatible.
c          = -14: DSAUPD  did not find any eigenvalues to sufficient
c                 accuracy.
c          = -15: HOWMNY must be one of 'A' or 'S' if RVEC = .true.
c          = -16: HOWMNY = 'S' not yet implemented
c          = -17: DSEUPD  got a different count of the number of converged
c                 Ritz values than DSAUPD  got.  This indicates the user
c                 probably made an error in passing data from DSAUPD  to
c                 DSEUPD  or that the data was modified before entering
c                 DSEUPD .

	*/
	int INFO;


	/*
	RVEC    LOGICAL  (INPUT) 
c          Specifies whether Ritz vectors corresponding to the Ritz value 
c          approximations to the eigenproblem A*z = lambda*B*z are computed.
c
c             RVEC = .FALSE.     Compute Ritz values only.
c
c             RVEC = .TRUE.      Compute Ritz vectors.
	*/
	bool RVEC;

	/*
	HOWMNY  Character*1  (INPUT) 
c          Specifies how many Ritz vectors are wanted and the form of Z
c          the matrix of Ritz vectors. See remark 1 below.
c          = 'A': compute NEV Ritz vectors;
c          = 'S': compute some of the Ritz vectors, specified
c                 by the logical array SELECT.
	*/
	char HOWMNY;



	/*
	SELECT  Logical array of dimension NCV.  (INPUT/WORKSPACE)
c          If HOWMNY = 'S', SELECT specifies the Ritz vectors to be
c          computed. To select the Ritz vector corresponding to a
c          Ritz value D(j), SELECT(j) must be set to .TRUE.. 
c          If HOWMNY = 'A' , SELECT is used as a workspace for
c          reordering the Ritz values.
	*/
	int *SELECT;


	/*
	D       Double precision  array of dimension NEV.  (OUTPUT)
c          On exit, D contains the Ritz value approximations to the
c          eigenvalues of A*z = lambda*B*z. The values are returned
c          in ascending order. If IPARAM(7) = 3,4,5 then D represents
c          the Ritz values of OP computed by dsaupd  transformed to
c          those of the original eigensystem A*z = lambda*B*z. If 
c          IPARAM(7) = 1,2 then the Ritz values of OP are the same 
c          as the those of A*z = lambda*B*z.
	*/
	double *D;


	/*
	 Z       Double precision  N by NEV array if HOWMNY = 'A'.  (OUTPUT)
c          On exit, Z contains the B-orthonormal Ritz vectors of the
c          eigensystem A*z = lambda*B*z corresponding to the Ritz
c          value approximations.
c          If  RVEC = .FALSE. then Z is not referenced.
c          NOTE: The array Z may be set equal to first NEV columns of the 
c          Arnoldi/Lanczos basis array V computed by DSAUPD .
	*/
	double *Z;


	/*
	LDZ     Integer.  (INPUT)
c          The leading dimension of the array Z.  If Ritz vectors are
c          desired, then  LDZ .ge.  max( 1, N ).  In any case,  LDZ .ge. 1.
	*/
	int LDZ;

	/*
c  SIGMA   Double precision   (INPUT)
c          If IPARAM(7) = 3,4,5 represents the shift. Not referenced if
c          IPARAM(7) = 1 or 2.
	*/
	double SIGMA;



};

