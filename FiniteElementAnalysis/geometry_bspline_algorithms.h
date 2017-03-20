#ifndef _GEOMETRY_BSPLINE_ALGORITHMS_H_
#define _GEOMETRY_BSPLINE_ALGORITHMS_H_


namespace geometry {

	const double double_epsilon = DBL_EPSILON; //1e-16;



	/*Basis function value at x*/
	static double deboor_value(int i, int p, double *knot, int lknot, double x){

		if (i < 0) throw "deboor error: require i >= 0";
		if (p < 0) throw "deboor error: require p >= 0";
		if (i + p + 1 >= lknot) throw "deboor error: require i + p + 1 < knot.length";

		if (x == knot[lknot - p - 1]){ x -= x * double_epsilon; }

		if (p == 0){
			if (knot[i] <= x && x < knot[i + 1]){
				return 1;
			}
			else{
				return 0;
			}
		}
		double d1 = knot[i + p] - knot[i];
		double d2 = knot[i + p + 1] - knot[i + 1];

		double val1 = 0;
		double val2 = 0;
		if (d1 > 0){
			val1 = (x - knot[i]) / d1 * deboor_value(i, p - 1, knot, lknot, x);
		}
		if (d2 > 0) {
			val2 = (knot[i + p + 1] - x) / d2 * deboor_value(i + 1, p - 1, knot, lknot, x);
		}
		return val1 + val2;
	}


	/*Basis function derivative value at x*/
	static double deboor_derivative(int i, int p, double *knot, int lknot, double x){

		if (i < 0) throw "deboor error: require i >= 0";
		if (p < 1) throw "deboor error: require p >= 1";
		if (i + p + 1 >= lknot) throw "deboor error: require i + p + 1 < knot.length";

		double value1 = deboor_value(i, p - 1, knot, lknot, x);
		double value2 = deboor_value(i + 1, p - 1, knot, lknot, x);

		double d1 = knot[i + p] - knot[i];
		double d2 = knot[i + p + 1] - knot[i + 1];

		double val1 = 0;
		double val2 = 0;

		if (d1 > 0){
			val1 += p / d1 * value1;
		}
		if (d2 > 0) {
			val2 += -p / d2 * value2;
		}
		return val1 + val2;
	}


	/*Basis function value and derivative value at x*/
	static void deboor_values(int i, int p, double *knot, int lknot, double x, double *result){

		double value1 = deboor_value(i, p - 1, knot, lknot, x);
		double value2 = deboor_value(i + 1, p - 1, knot, lknot, x);

		double d1 = knot[i + p] - knot[i];
		double d2 = knot[i + p + 1] - knot[i + 1];

		for (int j = 0; j < p; j++) result[j] = 0;
		if (d1 > 0){
			result[0] += (x - knot[i]) / d1 * value1;
			result[1] += p / d1 * value1;
		}
		if (d2 > 0) {
			result[0] += (knot[i + p + 1] - x) / d2 * value2;
			result[1] += -p / d2 * value2;
		}

	}




	/*
	Evaluates a bspline curve at parameter t using deboor algorithm
	*/
	template <class Type>
	static Type evaluate_curve(int p, const double* knot, int lknot, const Type* polygon, double t){
		if (t == knot[lknot - p - 1]){ t -= t * double_epsilon; }

		int k = find_span(lknot, knot, t);
		return evaluate_curve_deboor(p, k, p, knot, lknot, polygon, t);
	}


	template <class Type>
	static Type evaluate_curve_deboor(int p, int i, int j, const double *knot, int lknot, const Type *polygon, double t){
		if (j == 0) return polygon[i];
		double d = (knot[i + p + 1 - j] - knot[i]);
		if (d > 0)
		{
			double a = (t - knot[i]) / d;
			Type res1 = evaluate_curve_deboor(p, i - 1, j - 1, knot, lknot, polygon, t);
			res1 *= (1 - a);
			Type res2 = evaluate_curve_deboor(p, i, j - 1, knot, lknot, polygon, t);
			res2 *= a;
			res1 += res2;
			return res1;
			//return evaluate_curve_deboor(p, i - 1, j - 1, knot, lknot, polygon, t) * (1 - a) + evaluate_curve_deboor(p, i, j - 1, knot, lknot, polygon, t) * a;
		}
	}







	/*
	Oslo algorithm taken from Cohen et al.
	*/
	static int find_span(int KN, const double *TAU, double Tj){
		for (int i = 0; i < KN; i++){
			if (TAU[i] <= Tj && Tj < TAU[i + 1]) return i;
		}
	}

	template <class Type>
	static Type olso_subdiv(Type *P, int K, double *TAU, double *T, int RPI, int i, int j){
		int r = RPI - 1;
		if (r > 0){
			Type PP2; PP2 = 0.0;
			Type PP1; PP1 = 0.0;
			double P1 = T[j + K - r] - TAU[i];
			double P2 = TAU[i + K - r] - T[j + K - r];
			if (P1 != 0) PP1 = olso_subdiv(P, K, TAU, T, r, i, j);
			if (P2 != 0) PP2 = olso_subdiv(P, K, TAU, T, r, i - 1, j);
			if (abs(P1 + P2) > 1e-15){
				PP1 *= P1 / (P1 + P2);
				PP2 *= P2 / (P1 + P2);
				PP1 += PP2;
				return PP1;
				//return (P1 * PP1 + P2 * PP2) / (P1 + P2);
			}
			else{
				return Type();
			}
		}
		else{
			return P[i];
		}
	}

	/*
	N = where original polygon has N+1 vertices
	P = original polygon of Type
	K = order of bspline curve (i.e. p+1)
	TAU = original knot vector of length N + K
	T = refinement knot vector of length Q > N + K
	D = subdivided polygon of Type
	*/
	template <class Type>
	static void olso_insertion(int N, Type *P, int K, double* TAU, double *T, int Q, Type *D){
		int MU;
		for (int j = 0; j <= Q - K; j++){
			MU = find_span(K + N, TAU, T[j]);
			D[j] = olso_subdiv<Type>(P, K, TAU, T, K, MU, j);
		}
	}







}


#endif _GEOMETRY_BSPLINE_ALGORITHMS_H_