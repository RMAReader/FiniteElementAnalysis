#ifndef _OPTIMIZER_H_
#define _OPTIMIZER_H_

//#ifdef FEADLL_EXPORTS
//#define FEADLL_API __declspec(dllexport) 
//#else
//#define FEADLL_API __declspec(dllimport) 
//#endif

#include "dlib\optimization.h"
#include "BSplineSolid.h"
#include "arpack_shiftinvert.h"
#include "Element_BSplineSolid.h"
#include <iostream>


using namespace std;
using namespace dlib;


// ----------------------------------------------------------------------------------------

// In dlib, the general purpose solvers optimize functions that take a column
// vector as input and return a double.  So here we make a typedef for a
// variable length column vector of doubles.  This is the type we will use to
// represent the input to our objective functions which we will be minimizing.
typedef matrix<double, 0, 1> column_vector;

// ----------------------------------------------------------------------------------------
// Below we create a few functions.  When you get down into main() you will see that
// we can use the optimization algorithms to find the minimums of these functions.
// ----------------------------------------------------------------------------------------

//double rosen(const column_vector& m)
///*
//This function computes what is known as Rosenbrock's function.  It is
//a function of two input variables and has a global minimum at (1,1).
//So when we use this function to test out the optimization algorithms
//we will see that the minimum found is indeed at the point (1,1).
//*/
//{
//	const double x = m(0);
//	const double y = m(1);
//
//	// compute Rosenbrock's function and return the result
//	return 100.0*pow(y - x*x, 2) + pow(1 - x, 2);
//}
//
//// This is a helper function used while optimizing the rosen() function.  
//const column_vector rosen_derivative(const column_vector& m)
///*!
//ensures
//- returns the gradient vector for the rosen function
//!*/
//{
//	const double x = m(0);
//	const double y = m(1);
//
//	// make us a column vector of length 2
//	column_vector res(2);
//
//	// now compute the gradient vector
//	res(0) = -400 * x*(y - x*x) - 2 * (1 - x); // derivative of rosen() with respect to x
//	res(1) = 200 * (y - x*x);              // derivative of rosen() with respect to y
//	return res;
//}
//
//// This function computes the Hessian matrix for the rosen() fuction.  This is
//// the matrix of second derivatives.
//matrix<double> rosen_hessian(const column_vector& m)
//{
//	const double x = m(0);
//	const double y = m(1);
//
//	matrix<double> res(2, 2);
//
//	// now compute the second derivatives 
//	res(0, 0) = 1200 * x*x - 400 * y + 2; // second derivative with respect to x
//	res(1, 0) = res(0, 1) = -400 * x;   // derivative with respect to x and y
//	res(1, 1) = 200;                 // second derivative with respect to y
//	return res;
//}
//
//// ----------------------------------------------------------------------------------------
//
//class test_function
//{
//	/*
//	This object is an example of what is known as a "function object" in C++.
//	It is simply an object with an overloaded operator().  This means it can
//	be used in a way that is similar to a normal C function.  The interesting
//	thing about this sort of function is that it can have state.
//
//	In this example, our test_function object contains a column_vector
//	as its state and it computes the mean squared error between this
//	stored column_vector and the arguments to its operator() function.
//
//	This is a very simple function, however, in general you could compute
//	any function you wanted here.  An example of a typical use would be
//	to find the parameters of some regression function that minimized
//	the mean squared error on a set of data.  In this case the arguments
//	to the operator() function would be the parameters of your regression
//	function.  You would loop over all your data samples and compute the output
//	of the regression function for each data sample given the parameters and
//	return a measure of the total error.   The dlib optimization functions
//	could then be used to find the parameters that minimized the error.
//	*/
//public:
//
//	test_function(
//		const column_vector& input
//		)
//	{
//		target = input;
//	}
//
//	double operator() (const column_vector& arg) const
//	{
//		// return the mean squared error between the target vector and the input vector
//		return mean(squared(target - arg));
//	}
//
//private:
//	column_vector target;
//};
//
//// ----------------------------------------------------------------------------------------
//
//class rosen_model
//{
//	/*!
//	This object is a "function model" which can be used with the
//	find_min_trust_region() routine.
//	!*/
//
//public:
//	typedef ::column_vector column_vector;
//	typedef matrix<double> general_matrix;
//
//	double operator() (
//		const column_vector& x
//		) const {
//		return rosen(x);
//	}
//
//	void get_derivative_and_hessian(
//		const column_vector& x,
//		column_vector& der,
//		general_matrix& hess
//		) const
//	{
//		der = rosen_derivative(x);
//		hess = rosen_hessian(x);
//	}
//};
//
//// ----------------------------------------------------------------------------------------



class Optimizer
{
public:
	Optimizer();
	~Optimizer();
	void calibrate(BSplineSolid* in_solid, double in_density, double* in_D, int*, double* measured_modes, int number_modes);
	void tune_plate_harris(BSplineSolid* basis_solids, int basis_dimension, double in_density, double* in_D, double *target);

	int test();

};

class target_calibrate_elasticity
{
public:
	target_calibrate_elasticity(BSplineSolid* , double , double*, int*, double*, int , bool);
	~target_calibrate_elasticity();
	double operator()(const column_vector& arg) const;

private:
	BSplineSolid* solid;
	double density;
	double* D;
	int *target_modes;
	double *target_frequencies;
	int N;
	bool verbose;
};


class target_tune_belly_Harris
{
public:
	target_tune_belly_Harris(BSplineSolid*, int, double, double*,  double*, bool);
	~target_tune_belly_Harris();
	double operator()(const column_vector& arg) const;

private:
	BSplineSolid *basis_solids;
	int basis_dimension;
	double density;
	double* D;
	double *target;
	bool verbose;
};


#endif _OPTIMIZER_H_