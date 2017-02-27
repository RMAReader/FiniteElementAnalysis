#include "stdafx.h"
#include "Optimizer.h"



Optimizer::Optimizer()
{
}


Optimizer::~Optimizer()
{
}

void Optimizer::calibrate(BSplineSolid* in_solid, double in_density, double* in_D, int *measured_modes, double *measured_freq, int number_modes){

	column_vector starting_point(9);
	column_vector lower_bound(9);
	column_vector upper_bound(9);

	//set starting point to vector of 9 elasticities, initial scaling of the input elasticities, and initial lagrange multiplier values
	for (int i = 0; i < 9; i++){
		starting_point(i) = 0;	//elasticities scaling
		lower_bound(i) = -1e100;
		upper_bound(i) = 1e100;
	}
	
	for (int i = 0; i < 9; i++){
		cout << "starting_point(" << i << ") = " << starting_point(i) << endl;
	}

	try{
		find_min_bobyqa(target_calibrate_elasticity(in_solid, in_density, in_D, measured_modes, measured_freq, number_modes, true),
			starting_point,
			20,    // number of interpolation points
			lower_bound,  // lower bound constraint
			upper_bound,  // upper bound constraint
			2.0,    // initial trust region radius
			1e-6,  // stopping trust region radius
			100    // max number of objective function evaluations
			);
		cout << "test_function solution:\n" << starting_point << endl;
		//copy result into elasticity array D
		for (int i = 0; i < 9; i++){
			in_D[i] *= exp(starting_point(i));
		}

	}
	catch (bobyqa_failure e){
		cout << e.what() << endl;
	}
}

void Optimizer::tune_plate_harris(BSplineSolid* basis_solids, int basis_dimension,double in_density, double* in_D, double *target){

	column_vector starting_point(basis_dimension - 1);
	column_vector lower_bound(basis_dimension - 1);
	column_vector upper_bound(basis_dimension - 1);

	//set starting point to vector of 1 scaling value
	for (int i = 0; i < basis_dimension - 1; i++){
		starting_point(i) = 0.5;	
		lower_bound(i) = 0;
		upper_bound(i) = 1;
	}

	for (int i = 0; i < basis_dimension - 1; i++){
		cout << "starting_point(" << i << ") = " << starting_point(i) << endl;
	}

	try{
		find_min_bobyqa(target_tune_belly_Harris(basis_solids, basis_dimension, in_density, in_D, target, true),
			starting_point,
			2 * basis_dimension,    // number of interpolation points
			lower_bound,  // lower bound constraint
			upper_bound,  // upper bound constraint
			0.3,    // initial trust region radius
			1e-6,  // stopping trust region radius
			100    // max number of objective function evaluations
			);
		cout << "test_function solution:\n" << starting_point << endl;
		
	}
	catch (bobyqa_failure e){
		cout << e.what() << endl;
	}
}


target_calibrate_elasticity::target_calibrate_elasticity(BSplineSolid* in_solid, double in_density, double* in_D, int* modes, double* measured_freq, int number_modes, bool in_verbose)
{
	solid = in_solid;
	density = in_density;
	D = in_D;
	target_modes = modes;
	target_frequencies = measured_freq;
	N = number_modes;
	verbose = in_verbose;
}
target_calibrate_elasticity::~target_calibrate_elasticity()
{
}
double target_calibrate_elasticity::operator()(const column_vector& arg) const
{

	for (int i = 0; i < 9; i++){
		cout << "arg(" << i << ") = " << arg(i) << endl;
	}

	double newD[9];
	for (int i = 0; i < 9; i++){
		newD[i] = exp(arg(i)) * D[i];
	}

	Element_BSplineSolid element(*solid, newD, density);
	arpack_shiftinvert solver(element.getM(), element.getK());
	solver.shift(0);
	solver.factorizeLU();
	solver.solve();

	//squared log-ratio of simulated frequencies to measured frequencies
	double value = 0;
	for (int i = 0; i < N; i++){
		double inc = log(target_frequencies[i] / solver.get_frequency(target_modes[i]));
		value += inc * inc;
	}

	if (verbose){
		ofstream log("C:\\Private\\calibrationLog.csv", ios::out | ios::app);

		if (log.is_open()){
			time_t rawtime;
			struct tm * timeinfo;
			char buffer[80];
			time(&rawtime);
			timeinfo = localtime(&rawtime);
			strftime(buffer, 80, "%c", timeinfo);

			log << buffer << ",";
			log << value << ",";
			for (int i = 0; i < 9; i++){
				log << arg(i) << ",";
			}
			for (int i = 0; i < 9; i++){
				log << exp(arg(i)) * D[i] << ",";
			}
			for (int i = 0; i < 20; i++){
				log << solver.get_frequency(i) << ",";
			}
			log << endl;
		}
		else{
			cout << "failed to write to calibration log" << endl;
		}
	}
	return value;
}




target_tune_belly_Harris::target_tune_belly_Harris(BSplineSolid* basis_solids, int basis_dimension, double density, double* D, double* target, bool verbose)
{
	this->basis_solids = basis_solids;
	this->basis_dimension = basis_dimension;
	this->density = density;
	this->D = D;
	this->target = target;
	this->verbose = verbose;
}
target_tune_belly_Harris::~target_tune_belly_Harris()
{
}
/*
Harris method of tuning Belly:
Optimise plate thicknesses so that K and Ratio equal target values, where
 - stiffness factor K = mass * 0.25 * (mode2 + mode5)^2
 - Ratio  = mode5/mode2
Thickness of plate is varied uniformly, i.e. taking constant proportion of thickness basis 1
Harris suggests K = 4250 and Ratio = 2
*/
double target_tune_belly_Harris::operator()(const column_vector& arg) const
{
	BSplineSolid newSolid;
	newSolid = basis_solids[0];
	for (int i = 0; i < newSolid.npoints(); i++){
		for (int j = 0; j < basis_dimension - 1; j++){
			newSolid.cz[i] += arg(j) * (basis_solids[j + 1].cz[i] - basis_solids[0].cz[i]);
		}
	}

	Element_BSplineSolid element(newSolid, D, density);
	arpack_shiftinvert solver(element.getM(), element.getK());

	solver.shift(0);
	solver.factorizeLU();
	solver.solve();

	//squared log-ratio of target values to simulated values
	double f[2];
	f[0] = element.get_mass() * 0.25 * (solver.get_frequency(1) + solver.get_frequency(4))* (solver.get_frequency(1) + solver.get_frequency(4));
	f[1] = solver.get_frequency(4) / solver.get_frequency(1);
	
	double value = 0;
	for (int i = 0; i < 2; i++){
		double inc = log(target[i] / f[i]);
		value += inc * inc;
	}

	if (verbose){
		ofstream log("C:\\Private\\tuningLog.csv", ios::out | ios::app);

		if (log.is_open()){
			time_t rawtime;
			struct tm * timeinfo;
			char buffer[80];
			time(&rawtime);
			timeinfo = localtime(&rawtime);
			strftime(buffer, 80, "%c", timeinfo);

			log << buffer << ",";
			log << value << ",";
			for (int i = 0; i < 1; i++){
				log << arg(i) << ",";
			}
			for (int i = 0; i < 2; i++){
				log << f[i] << ",";
			}
			for (int i = 0; i < 20; i++){
				log << solver.get_frequency(i) << ",";
			}
			log << endl;
		}
		else{
			cout << "failed to write to calibration log" << endl;
		}
	}
	return value;
}



//
//
//
//int Optimizer::test(){
//	try
//	{
//		// make a column vector of length 2
//		column_vector starting_point(2);
//
//
//		// Set the starting point to (4,8).  This is the point the optimization algorithm
//		// will start out from and it will move it closer and closer to the function's 
//		// minimum point.   So generally you want to try and compute a good guess that is
//		// somewhat near the actual optimum value.
//		starting_point = 4, 8;
//
//		// The first example below finds the minimum of the rosen() function and uses the
//		// analytical derivative computed by rosen_derivative().  Since it is very easy to
//		// make a mistake while coding a function like rosen_derivative() it is a good idea
//		// to compare your derivative function against a numerical approximation and see if
//		// the results are similar.  If they are very different then you probably made a 
//		// mistake.  So the first thing we do is compare the results at a test point: 
//		cout << "Difference between analytic derivative and numerical approximation of derivative: "
//			<< length(derivative(rosen)(starting_point)-rosen_derivative(starting_point)) << endl;
//
//
//		cout << "Find the minimum of the rosen function()" << endl;
//		// Now we use the find_min() function to find the minimum point.  The first argument
//		// to this routine is the search strategy we want to use.  The second argument is the 
//		// stopping strategy.  Below I'm using the objective_delta_stop_strategy which just 
//		// says that the search should stop when the change in the function being optimized 
//		// is small enough.
//
//		// The other arguments to find_min() are the function to be minimized, its derivative, 
//		// then the starting point, and the last is an acceptable minimum value of the rosen() 
//		// function.  That is, if the algorithm finds any inputs to rosen() that gives an output 
//		// value <= -1 then it will stop immediately.  Usually you supply a number smaller than 
//		// the actual global minimum.  So since the smallest output of the rosen function is 0 
//		// we just put -1 here which effectively causes this last argument to be disregarded.
//
//		find_min(bfgs_search_strategy(),  // Use BFGS search algorithm
//			objective_delta_stop_strategy(1e-7), // Stop when the change in rosen() is less than 1e-7
//			rosen, rosen_derivative, starting_point, -1);
//		// Once the function ends the starting_point vector will contain the optimum point 
//		// of (1,1).
//		cout << "rosen solution:\n" << starting_point << endl;
//
//
//		// Now let's try doing it again with a different starting point and the version
//		// of find_min() that doesn't require you to supply a derivative function.  
//		// This version will compute a numerical approximation of the derivative since 
//		// we didn't supply one to it.
//		starting_point = -94, 5.2;
//		find_min_using_approximate_derivatives(bfgs_search_strategy(),
//			objective_delta_stop_strategy(1e-7),
//			rosen, starting_point, -1);
//		// Again the correct minimum point is found and stored in starting_point
//		cout << "rosen solution:\n" << starting_point << endl;
//
//
//		// Here we repeat the same thing as above but this time using the L-BFGS 
//		// algorithm.  L-BFGS is very similar to the BFGS algorithm, however, BFGS 
//		// uses O(N^2) memory where N is the size of the starting_point vector.  
//		// The L-BFGS algorithm however uses only O(N) memory.  So if you have a 
//		// function of a huge number of variables the L-BFGS algorithm is probably 
//		// a better choice.
//		starting_point = 0.8, 1.3;
//		find_min(lbfgs_search_strategy(10),  // The 10 here is basically a measure of how much memory L-BFGS will use.
//			objective_delta_stop_strategy(1e-7).be_verbose(),  // Adding be_verbose() causes a message to be 
//			// printed for each iteration of optimization.
//			rosen, rosen_derivative, starting_point, -1);
//
//		cout << endl << "rosen solution: \n" << starting_point << endl;
//
//		starting_point = -94, 5.2;
//		find_min_using_approximate_derivatives(lbfgs_search_strategy(10),
//			objective_delta_stop_strategy(1e-7),
//			rosen, starting_point, -1);
//		cout << "rosen solution: \n" << starting_point << endl;
//
//
//
//
//		// dlib also supports solving functions subject to bounds constraints on
//		// the variables.  So for example, if you wanted to find the minimizer
//		// of the rosen function where both input variables were in the range
//		// 0.1 to 0.8 you would do it like this:
//		starting_point = 0.1, 0.1; // Start with a valid point inside the constraint box.
//		find_min_box_constrained(lbfgs_search_strategy(10),
//			objective_delta_stop_strategy(1e-9),
//			rosen, rosen_derivative, starting_point, 0.1, 0.8);
//		// Here we put the same [0.1 0.8] range constraint on each variable, however, you
//		// can put different bounds on each variable by passing in column vectors of
//		// constraints for the last two arguments rather than scalars.  
//
//		cout << endl << "constrained rosen solution: \n" << starting_point << endl;
//
//		// You can also use an approximate derivative like so:
//		starting_point = 0.1, 0.1;
//		find_min_box_constrained(bfgs_search_strategy(),
//			objective_delta_stop_strategy(1e-9),
//			rosen, derivative(rosen), starting_point, 0.1, 0.8);
//		cout << endl << "constrained rosen solution: \n" << starting_point << endl;
//
//
//
//
//		// In many cases, it is useful if we also provide second derivative information
//		// to the optimizers.  Two examples of how we can do that are shown below.  
//		starting_point = 0.8, 1.3;
//		find_min(newton_search_strategy(rosen_hessian),
//			objective_delta_stop_strategy(1e-7),
//			rosen,
//			rosen_derivative,
//			starting_point,
//			-1);
//		cout << "rosen solution: \n" << starting_point << endl;
//
//		// We can also use find_min_trust_region(), which is also a method which uses
//		// second derivatives.  For some kinds of non-convex function it may be more
//		// reliable than using a newton_search_strategy with find_min().
//		starting_point = 0.8, 1.3;
//		find_min_trust_region(objective_delta_stop_strategy(1e-7),
//			rosen_model(),
//			starting_point,
//			10 // initial trust region radius
//			);
//		cout << "rosen solution: \n" << starting_point << endl;
//
//
//
//
//		// Now let's look at using the test_function object with the optimization 
//		// functions.  
//		cout << "\nFind the minimum of the test_function" << endl;
//
//		column_vector target(4);
//		starting_point.set_size(4);
//
//		// This variable will be used as the target of the test_function.   So,
//		// our simple test_function object will have a global minimum at the
//		// point given by the target.  We will then use the optimization 
//		// routines to find this minimum value.
//		target = 3, 5, 1, 7;
//
//		// set the starting point far from the global minimum
//		starting_point = 1, 2, 3, 4;
//		find_min_using_approximate_derivatives(bfgs_search_strategy(),
//			objective_delta_stop_strategy(1e-7),
//			test_function(target), starting_point, -1);
//		// At this point the correct value of (3,5,1,7) should be found and stored in starting_point
//		cout << "test_function solution:\n" << starting_point << endl;
//
//		// Now let's try it again with the conjugate gradient algorithm.
//		starting_point = -4, 5, 99, 3;
//		find_min_using_approximate_derivatives(cg_search_strategy(),
//			objective_delta_stop_strategy(1e-7),
//			test_function(target), starting_point, -1);
//		cout << "test_function solution:\n" << starting_point << endl;
//
//
//
//		// Finally, let's try the BOBYQA algorithm.  This is a technique specially
//		// designed to minimize a function in the absence of derivative information.  
//		// Generally speaking, it is the method of choice if derivatives are not available.
//		starting_point = -4, 5, 99, 3;
//		find_min_bobyqa(test_function(target),
//			starting_point,
//			9,    // number of interpolation points
//			uniform_matrix<double>(4, 1, -1e100),  // lower bound constraint
//			uniform_matrix<double>(4, 1, 1e100),   // upper bound constraint
//			10,    // initial trust region radius
//			1e-6,  // stopping trust region radius
//			100    // max number of objective function evaluations
//			);
//		cout << "test_function solution:\n" << starting_point << endl;
//		return 0;
//	}
//	catch (std::exception& e)
//	{
//		cout << e.what() << endl;
//		return 1;
//	}
//}