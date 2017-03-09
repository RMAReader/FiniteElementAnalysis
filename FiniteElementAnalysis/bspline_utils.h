
#ifndef _BSPLINE_H_
#define _BSPLINE_H_

#include "iostream"
#include <exception>
#include <memory>
#include <vector>


#define VEC2F bspline::vec2<float>
#define VEC3F bspline::vec3<float>
#define VEC2D bspline::vec2<double>
#define VEC3D bspline::vec3<double>
#define CURVE2F bspline::curve<bspline::vec2<float>>
#define CURVE3F bspline::curve<bspline::vec3<float>>
#define CURVE2D bspline::curve<bspline::vec2<double>>
#define CURVE3D bspline::curve<bspline::vec3<double>>
#define LATTICE3F bspline::lattice<bspline::vec3<float>>
#define SURFACE3F bspline::surface<bspline::vec3<float>>
#define SURFACE3D bspline::surface<bspline::vec3<double>>
#define TRIANGLE3F bspline::triangle3f


namespace bspline{

	const double double_epsilon = DBL_EPSILON; //1e-16;
	const double double_min_value = DBL_MIN;
	const double double_max_value = DBL_MAX;

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
	static Type evaluate_curve(int p, const double &knot, int lknot, const Type &polygon, double t){
		if (t == (&knot)[lknot - p - 1]){ t -= t * double_epsilon; }

		int k = find_span(lknot, &knot, t);
		return evaluate_curve_deboor(p, k, p, &knot, lknot, &polygon, t);
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
			Type PP2; PP2 = 0;
			Type PP1; PP1 = 0;
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


	/*Returns the number of knot spans covered by a knot vector*/
	static int number_spans(int p, int lknot){
		return lknot - 2 * p - 1;
	}


	/*creates a new knot vector with each span split evenly into n spans*/
	static void split_knots(int p, double* input_knot, int lknot, int n, double* output_knot, int output_lknot){
		if (n < 1) return;
		int spans = number_spans(p, lknot);
		if (spans < 1)return;

		/*output_lknot = lknot + spans * (n - 1);
		output_knot = new double[output_lknot];*/

		for (int i = 0; i <= p; i++){
			output_knot[i] = input_knot[i];
			output_knot[output_lknot - i - 1] = input_knot[lknot - i - 1];
		}
		for (int i = 0; i < spans; i++){
			for (int j = 0; j < n; j++){
				output_knot[p + n * i + j] = (double)j / n * input_knot[i + 1] + (1 - (double)j / n)*input_knot[i];
			}
		}

	}

	static void merge_knots(double* knot, int m, double* extra_knot, int n, double** new_knot, int& k)
	{
		double* tmp_knot = new double[n + m];
		int j = 0;
		k = 0;
		
		for (int i = 0; i < m; i++)
		{
			tmp_knot[k] = knot[i];
			k++;
			if (i < m - 1){
				while (extra_knot[j] < knot[i]) j++;
				while (extra_knot[j] < knot[i + 1]){
					tmp_knot[k] = extra_knot[j];
					j++;
					k++;
				}
			}
		}
		k--;
		*new_knot = new double[k];
		for (int i = 0; i < k; i++){
			(*new_knot)[i] = tmp_knot[i];
		}
		delete[] tmp_knot;
	}
	









	template <class Type>
	class curve {
		int _p;
		int _n;
		double *knot;
		int _lknot;
		Type *points;

	public:
		curve()
		{
		}

		curve(int p, int n, int lknot, double* knot, Type* points)
		{
			_p = p;
			_n = n;
			_lknot = lknot;
			this->knot = new double[lknot];
			for (int i = 0; i < lknot; i++)
			{
				this->knot[i] = knot[i];
			}
			this->points = new Type[n];
			for (int i = 0; i < n; i++)
			{
				this->points[i] = points[i];
			}
		}

		curve(int order, Type point)
		{
			_p = order;
			_n = 1;
			_lknot = _p + _n + 1;
			points = new Type[_n];
			points[0] = point;
			knot = new double[_lknot];
			for (int i = 0; i < _lknot; i++){
				knot[i] = fmin(fmax((double)(i - _p) / (_n - _p), 0), 1);
			}
		}
		
		curve(const curve<Type> &c)
		{
			_p = c._p;
			_n = c._n;
			_lknot = c._lknot;
			knot = new double[_lknot];
			for (int i = 0; i < _lknot; i++){
				knot[i] = c.knot[i];
			}
			points = new Type[_n];
			for (int i = 0; i < _n; i++){
				points[i] = c.points[i];
			}
		}

		virtual ~curve()
		{
			delete[] knot;
			delete[] points;
		}

		void initialise(int p, int n, int lknot, double* knot, Type* points)
		{
			_p = p;
			_n = n;
			_lknot = lknot;
			this->knot = knot;
			this->points = points;
		}

		Type evaluate(double t)
		{
			return evaluate_curve(_p, *knot, _lknot, *points, t);
		}

		void append(Type point)
		{
			_n += 1;
			_lknot += 1;

			Type *newPoints = new Type[_n];
			for (int i = 0; i < _n - 1; i++){
				newPoints[i] = points[i];
			}
			newPoints[_n - 1] = point;
			delete[] points;
			points = newPoints;

			double *newknot = new double[_lknot];
			for (int i = 0; i < _lknot; i++){
				newknot[i] = fmin(fmax((double)(i - _p) / (_n - _p), 0),1);
			}
			delete[] knot;
			knot = newknot;
		}
		
		curve<Type>* insert_knots(double* extra_knot, int n)
		{
			//count extra_knots within range of current knot vector
			int m = 0;
			for (int i = 0; i < n; i++){
				if (extra_knot[i] >= knot[0] && extra_knot[i] <= knot[_lknot - 1]){ m++; }
			}
			double* new_knot = new double[m + _lknot];
			int last_j = 0;
			int k = 0;
			for (int i = 0; i < _lknot - 1; i++){
				new_knot[k] = knot[i];
				k++;
				for (int j = last_j; j < n; j++){
					if (extra_knot[j] >= knot[i] && extra_knot[j] <= knot[i + 1]){
						new_knot[k] = extra_knot[j];
						last_j++;
						k++;
					}
					else {
						if (extra_knot[j] < knot[i + 1]){ last_j++; }
						break;
					}
				}
			}
			new_knot[k] = knot[_lknot-1];
			Type* new_points = new Type[m + _n];

			bspline::olso_insertion(_n-1, points, _p+1, knot, new_knot,  m + _n + _p, new_points);
			curve<Type>* output = new curve<Type>(_p, _n + m, _lknot + m, new_knot, new_points);

			delete[] new_points;
			delete[] new_knot;
			return output;
		}


		curve<Type>* trim_curve(double min_t, double max_t)
		{
			double* extra_knot=  new double[(_p+1) * 2];
			for (int i = 0; i <= _p; i++){
				extra_knot[i] = min_t;
				extra_knot[i + _p + 1] = max_t;
			}
			curve<Type>* c1 = this->insert_knots(extra_knot, (_p + 1) * 2);

			int i = 0;
			int j = c1->lKnot() - 1;
			while (c1->knot[i] <= min_t && i<c1->lKnot()-1) i++;
			while (c1->knot[j] >= max_t && j > 0) j--;

			i -= (_p + 1);
			j += (_p + 1);

			int new_lknot = j - i + 1; 
			double* new_knot = new double[new_lknot];
			
			int new_n = new_lknot - _p - 1;
			Type* new_points = new Type[new_n];

			for (int k = 0; k < new_lknot; k++){
				new_knot[k] = c1->getKnot(k + i);
			}
			for (int k = 0; k < new_n; k++){
				new_points[k] = c1->get(k + i);
			}
			return new curve<Type>(_p, new_n, new_lknot, new_knot, new_points);
		}

		float length(double minParam, double maxParam)
		{
			double min_t, max_t;
			if (maxParam > minParam){
				min_t = (this->minParam() < minParam) ? minParam : this->minParam();
				max_t = (this->maxParam() > maxParam) ? maxParam : this->maxParam();
			}
			else{
				//min_t = (this->maxParam() > minParam) ? minParam : this->maxParam();
				//max_t = (this->minParam() < maxParam) ? maxParam : this->minParam();
				max_t = (this->minParam() < minParam) ? minParam : this->minParam();
				min_t = (this->maxParam() > maxParam) ? maxParam : this->maxParam();
			}
			curve<Type>* c1 = this->trim_curve(min_t, max_t);
			
			float length1 = 0;
			for (int i = 0; i < c1->nPoints() - 1; i++){
				length1 += distance(c1->get(i), c1->get(i + 1));
			}

			delete c1;
			return length1;
		}


		void set(int i, Type point){
			if (i < 0 || i>_n - 1) throw "bspline curve error: require 0 <= i < _n";
			points[i] = point;
		}
		Type get(int i){
			if (i < 0 || i>_n - 1) throw "bspline curve error: require 0 <= i < _n";
			return points[i];
		}
		Type* item(int i){
			if (i < 0 || i>_n - 1) throw "bspline curve error: require 0 <= i < _n";
			return &points[i];
		}
		int nPoints(){
			return _n;
		}
		int getOrder(){
			return _p;
		}
		double getKnot(int i){
			return knot[i];
		}
		double* getKnotVector(){
			return knot;
		}
		int lKnot(){
			return _lknot;
		}
		double minParam(){
			return knot[_p];
		}
		double maxParam(){
			return knot[_lknot - _p - 1];
		}
	};


	

	template <class T>
	class lattice
	{
	public:
		std::vector<T> data;
		int rows, cols;

		lattice(){}

		lattice(lattice& base){
			rows = base.rows;
			cols = base.cols;
			data = base.data;
		}

		lattice(int rows, int cols){
			this->rows = rows;
			this->cols = cols;
			data = std::vector < T >(rows * cols);
			//points.resize(rows * cols);
		}
		
		lattice(std::vector<T>& p, int rows, int cols){
			if (p.size() == rows * cols){
				this->data = p;
				this->rows = rows;
				this->cols = cols;
			}
		}
		
		T& GetPoint(int row, int col){
			int i = row * cols + col;
			return data[i];
		}
		std::vector<T> GetRow(int row){
			std::vector < T > slice(cols);
			//slice.resize(cols);
			int offset = row * cols;
			for (int i = 0; i < cols; i++){
				slice[i] = data[offset + i];
			}
			return slice;
		}
		std::vector<T> GetCol(int col){
			std::vector < T > slice(rows);
			//slice.resize(cols);
			for (int i = 0; i < rows; i++){
				slice[i] = data[i * cols + col];
			}
			return slice;
		}
		void SetPoint(int row, int col, T p){
			data[row * cols + col] = p;
		}
	};


	template <class T>
	class surface
	{
	public:
		int p; 
		int q;
		lattice<T> points;
		std::vector<double> knotx;
		std::vector<double> knoty;

		surface(){}

		surface(surface& base){
			p = base.p;
			q = base.q;
			points = lattice<T>(base.points);
			knotx = base.knotx;
			knoty = base.knoty;
		}

		surface(int p, int q, std::vector<double>& knotx, std::vector<double>& knoty, lattice<T>& points)
		{
			this->p = p;
			this->q = q;
			this->points = points;
			this->knotx = knotx;
			this->knoty = knoty;
		}

		curve<T>* curve_x(double v){
			std::vector<T> c(points.rows);
			for (int j = 0; j < points.rows; j++){
				c[j] = evaluate_curve(q, *(knoty.data()), knoty.size(), *(points.GetCRow(j).data()), v);
			}
			return new curve<T>(p, c.size(), knotx.size(), knotx.data(), c.data());
		}

		curve<T>* curve_y(double u){
			std::vector<T> c(points.cols);
			for (int j = 0; j < points.cols; j++){
				c[j] = evaluate_curve(p, *(knotx.data()), knotx.size(), *(points.GetCol(j).data()), u);
			}
			return new curve<T>(q, c.size(), knoty.size(), knoty.data(), c.data());
		}

		T evaluate(double u, double v){
			std::vector<T> c(points.cols);
			for (int j = 0; j < points.cols; j++){
				c[j] = evaluate_curve(p, *(knotx.data()), knotx.size(), *(points.GetCol(j).data()), u);
			}
			return evaluate_curve(q, *(knoty.data()), knoty.size(), *(c.data()), v);
		}




		double minParam_x(){
			return knotx[p];
		}
		double maxParam_x(){
			return knotx[knotx.size() - p - 1];
		}
		double minParam_y(){
			return knoty[q];
		}
		double maxParam_y(){
			return knoty[knoty.size() - q - 1];
		}
	};




	template <class Type>
	class vec2
	{
	private:
		curve<Type> *parent;

	public:
		Type x;
		Type y;
	

		vec2()
		{
			x = 0;
			y = 0;
		}
		vec2(Type u, Type v)
		{
			x = u;
			y = v;
		}
		vec2(const vec2<Type> &v)
		{
			x = v.x;
			y = v.y;
		}

		virtual ~vec2(){}

		void set(Type x, Type y)
		{
			this->x = x;
			this->y = y;
		}
		void set(vec2& v)
		{
			this->x = v.x;
			this->y = v.y;
		}


		void operator=(vec2& v)
		{
			this->x = v.x;
			this->y = v.y;
		}
		void operator=(Type v)
		{
			this->x = v;
			this->y = v;
		}
		void operator +=(vec2& v)
		{
			this->x += v.x;
			this->y += v.y;
		}
		void operator -=(vec2& v)
		{
			this->x -= v.x;
			this->y -= v.y;
		}
		void operator *=(Type v)
		{
			this->x *= v;
			this->y *= v;
		}

	};
	
	template <class Type>
	class vec3
	{
	private:
		curve<Type> *parent;

	public:
		Type x;
		Type y;
		Type z;

		vec3()
		{
			x = 0;
			y = 0;
			z = 0;
		}
		vec3(Type u, Type v, Type w)
		{
			x = u;
			y = v;
			z = w;
		}
		vec3(const vec3<Type> &v)
		{
			x = v.x;
			y = v.y;
			z = v.z;
		}

		virtual ~vec3(){}

		void set(Type x, Type y, Type z)
		{
			this->x = x;
			this->y = y;
			this->z = z;
		}
		void set(vec3& v)
		{
			this->x = v.x;
			this->y = v.y;
			this->z = v.z;
		}


		void operator=(vec3& v)
		{
			this->x = v.x;
			this->y = v.y;
			this->z = v.z;
		}
		void operator=(Type v)
		{
			this->x = v;
			this->y = v;
			this->z = v;
		}
		void operator +=(vec3& v)
		{
			this->x += v.x;
			this->y += v.y;
			this->z += v.z;
		}
		void operator *=(Type v)
		{
			this->x *= v;
			this->y *= v;
			this->z *= v;
		}

	};


	//template <class Type>
	static float distance(vec2<float> v1, vec2<float> v2){
		vec2<float> d = v1;
		d -= v2;
		return sqrtf(d.x*d.x + d.y*d.y);
		return 0;
	}



	//template <class Type>
	static bool IntersectionParam(vec2<float>& p1, vec2<float>& p2, vec2<float>& q1, vec2<float>& q2, vec2<double>* param, bool verbose){
		float a = p1.x - p2.x;
		float b = p1.y - p2.y;
		float c = q1.x - q2.x;
		float d = q1.y - q2.y;

		double det = a * d - b * c;
		if (det == 0){ return false; }

		float e = q2.x - p2.x;
		float f = q2.y - p2.y;

		double alpha = (d * e - c * f) / det;
		double beta = (b * e - a * f) / det;

		param->x = alpha;
		param->y = beta;

		if (alpha >= 0 && alpha <= 1 && beta >= 0 && beta <= 1){

			if (verbose){
				std::cout << "intersection: " << std::endl;
				std::cout << "  p1=(" << p1.x << "," << p1.y << ") p2=(" << p2.x << "," << p2.y << ")" << std::endl;
				std::cout << "  q1=(" << q1.x << "," << q1.y << ") q2=(" << q2.x << "," << q2.y << ")" << std::endl;
				std::cout << "  alpha=" << alpha << ", beta=" << beta << std::endl;
			}
			return true;
		}
		else{
			return false;
		}
	}

	//returns true or false if c1 intersects c2.  If there is intersection, parameters of intersect are set.
	//template <class Type>
	static bool IntersectParam(curve<vec2<float>>& c1, curve<vec2<float>>& c2, double error, vec2<double>* param, bool verbose){
		
		int max_itr = 100;

		curve<vec2<float>> d1 = c1;
		curve<vec2<float>> d2 = c2;

		bool IntersectFound = false;
		bool RefineIntersect = true;
		int i, j;

		for (int itr = 1; itr <= max_itr; itr++){
			
			if (verbose) { std::cout << "iteration " << itr << std::endl; }

			for (i = 0; i < d1.nPoints() - 1; i++){
				for (j = 0; j < d2.nPoints() - 1; j++){
					IntersectFound = IntersectionParam(d1.get(i), d1.get(i + 1), d2.get(j), d2.get(j + 1), param, verbose);
					if (IntersectFound)	break;					
				}
				if (IntersectFound) break;
			}
			if (IntersectFound){
				double t1 = 0;
				for (int k = i; k < i + d1.getOrder()+1; k++){
					t1 += param->x * d1.getKnot(k) + (1 - param->x) * d1.getKnot(k + 1);
				}
				t1 /= (d1.getOrder() + 1);
				double t11[1] = { t1  };
				param->x = t1;

				double t2 = 0;
				for (int k = j; k < j + d2.getOrder()+1; k++){
					t2 += param->y * d2.getKnot(k) + (1 - param->y) * d2.getKnot(k + 1);
				}
				t2 /= (d2.getOrder() + 1);
				double t22[1] = { t2 };
				param->y = t2;

				vec2<float> di1 = d1.evaluate(param->x);
				vec2<float> di2 = d2.evaluate(param->y);

				d1 = *(d1.insert_knots(&(param->x), 1));
				d2 = *(d2.insert_knots(&(param->y), 1));

				if (verbose){
					std::cout << " d1(" << di1.x << "," << di1.y << ")";
					std::cout << " d2(" << di2.x << "," << di2.y << ")";
					std::cout << std::endl;
					std::cout << " distance = " << distance(di1, di2) << std::endl;
				}

				if (distance(di1, di2) < error) { return true; }

			}
			else{
				break;
			}
		}
		return false;
	}


	


	class offsetCurve2f
	{
	public:
		offsetCurve2f(curve<vec2<float>>*, double);
		void SetBase(curve<vec2<float>>*);
		curve<vec2<float>>* GetBase();
		void SetOffset(double x);
		double GetOffset();
		vec2<float> evaluate(double t);

	private:
		curve<vec2<float>>* base;
double offset;
	};



	class triangle3f
	{
	public:
		VEC3F p[3];

		triangle3f(){}

		triangle3f(VEC3F p0, VEC3F p1, VEC3F p2)
		{
			p[0] = p0;
			p[1] = p1;
			p[2] = p2;
		}

		float min_x()
		{
			float result = FLT_MAX;
			for (int i = 0; i < 3; i++){
				if (p[i].x < result){ result = p[i].x; }
			}
			return result;
		}
		float max_x()
		{
			float result = FLT_MIN;
			for (int i = 0; i < 3; i++){
				if (p[i].x > result){ result = p[i].x; }
			}
			return result;
		}
		float min_y()
		{
			float result = FLT_MAX;
			for (int i = 0; i < 3; i++){
				if (p[i].y < result){ result = p[i].y; }
			}
			return result;
		}
		float max_y()
		{
			float result = FLT_MIN;
			for (int i = 0; i < 3; i++){
				if (p[i].y > result){ result = p[i].y; }
			}
			return result;
		}
		float min_z()
		{
			float result = FLT_MAX;
			for (int i = 0; i < 3; i++){
				if (p[i].z < result){ result = p[i].z; }
			}
			return result;
		}
		float max_z()
		{
			float result = FLT_MIN;
			for (int i = 0; i < 3; i++){
				if (p[i].z > result){ result = p[i].z; }
			}
			return result;
		}
	};

	class triangle2f
	{
	public:
		VEC2F p[3];
	};




	class mesh3f
	{
	public:

		std::vector<VEC3F> points;
		std::vector<std::vector<int>> elements;
		
		lattice<std::vector<int>> nodes;
		std::vector<float> node_x;
		std::vector<float> node_y;

		inline triangle3f get_triangle(int index);
		VEC3F get_vertex(int element, int v);
		std::vector<triangle3f>* triangles();

		void get_mesh_range(VEC3F* min, VEC3F* max);
		inline void get_element_range(int i, VEC3F* min, VEC3F* max);
		void build_nodes(int nx, int ny);
	};




	static bool interior_circle_point(VEC2F& c, float r, VEC2F p1){
		p1 -= c;
		if (p1.x * p1.x + p1.y * p1.y < r * r){ return true; }
		return false;
	}


	//returns true if circle overlaps p1->p2 edge of triangle such that could overlap interior of triangle [p1, p2, p3]
	static bool interior_triangle_edge(VEC2F c, float r, VEC2F p1, VEC2F p2, VEC2F p3){

		//1. translate points so p1 is origin
		p2 -= p1;
		p3 -= p1;
		c -= p1;
		p1 -= p1;

		//2. rotate points about origin so p2 lies on positive x axis
		VEC2F q2(sqrtf(p2.x * p2.x + p2.y * p2.y), 0);
		
		VEC2F d(c.x * p3.x + c.y * p3.y, c.x * p3.y - c.y * p3.x);
		d *= 1 / q2.x;

		if (d.x + r < 0){ return false; }
		if (d.x - r > q2.x){ return false; }

		VEC2F q3(p2.x * p3.x + p2.y * p3.y, p2.x * p3.y - p2.y * p3.x);
		q3 *= 1 / q2.x;

		if (d.y + r < fminf(0, q3.y)){ return false; }
		if (d.y - r > fmaxf(0, q3.y)){ return false; }

		return true;
	}

	//returns true if circle and triangle overlap
	static bool interior_triangle(VEC2F& c, float r, VEC2F& p1, VEC2F& p2, VEC2F& p3){
		
		//1. check if any triangle vertices are inside cirlce (cheap)
		if (interior_circle_point(c, r, p1) == true){ return true; }
		if (interior_circle_point(c, r, p2) == true){ return true; }
		if (interior_circle_point(c, r, p3) == true){ return true; }

		//2. check if circle is on interior side of all triangle edges
		if (interior_triangle_edge(c, r, p1, p2, p3) == false){ return false; }
		if (interior_triangle_edge(c, r, p3, p1, p2) == false){ return false; }
		if (interior_triangle_edge(c, r, p2, p3, p1) == false){ return false; }

		return true;
	}

	
	static mesh3f* build_mesh(SURFACE3F* surface, int nx, int ny)
	{
		mesh3f* mesh = new mesh3f();

		CURVE3F* curve_y;

		float u, v;
		for (int i = 0; i < nx; i++){

			u = (1 - (float)i / (nx - 1)) * surface->minParam_x() + (float)i / (nx - 1) * surface->maxParam_x();
			curve_y = surface->curve_y(u);

			for (int j = 0; j < ny; j++){
				v = (1 - (float)j / (ny - 1)) * surface->minParam_y() + (float)j / (ny - 1) * surface->maxParam_y();

				mesh->points.push_back(curve_y->evaluate(v));

			}
		}

		int p1, p2, p3, p4, count;
		mesh->elements.reserve(2 * nx * ny);
		for (int i = ny; i < nx * ny; i+= ny){
			for (int j = 1; j < ny; j++){
				p1 = i + j;
				p2 = p1 - 1;
				p3 = p1 - ny;
				p4 = p3 - 1;
				mesh->elements.push_back({ p1, p2, p3 });
				mesh->elements.push_back({ p4, p3, p2 });
			}
		}

		return mesh;
	}


}


class base2f
{
public:
	base2f();
	virtual ~base2f();

private:
};



class point2f : public bspline::vec2 < float >
{

public:
	point2f();
	point2f(float u, float v);
	point2f(base2f* parent, float u, float v);

	void Translate(float u, float v);
	void Translate(point2f p);
	void SetParent(base2f* parent);
	base2f* GetParent();
	float GetDistance(point2f& point);

private:
	base2f *parent;

};

class point3f : public bspline::vec3 < float >
{
public:
	point3f();
	point3f(float u, float v, float w);
	point3f(const point3f &v);

};

class circle2f
{
public:
	circle2f();
	circle2f(point2f p, float d);
	circle2f(float x, float y, float d);
	circle2f(const circle2f &c);

	point2f centre;
	float diameter;
};



class curve2f : public base2f, public bspline::curve < point2f >
{
public:
	curve2f(int p, int n, int lknot, double* knot, point2f* points);
	curve2f(int order, point2f point);
	void append(point2f point);
	//curve2f* trim_curve(double min_t, double max_t);
};






#endif  _BSPLINE_H_