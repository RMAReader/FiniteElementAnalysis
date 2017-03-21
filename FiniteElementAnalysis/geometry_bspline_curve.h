#ifndef _GEOMETRY_BSPLINE_CURVE_H_
#define _GEOMETRY_BSPLINE_CURVE_H_


#include <vector>
//#include <algorithm>  
#include "geometry_bspline_knot.h"
#include "geometry_bspline_algorithms.h"
//#include "bspline_utils.h"

namespace geometry {
	namespace bspline{

		template <class T, class K>
		class curve {

		public:
			int _p = 0;
			geometry::bspline::knot<K> _knot;
			std::vector<T> _points;

			curve()	{}
			virtual ~curve(){}


			//copy assignment
			curve& operator=(const curve& other)
			{
				if (this != &other)
				{
					_p = other._p;
					_knot = other._knot;
					_points = other._points;
				}
				return *this;
			}
			//copy constructor
			curve(const curve& other)
			{
				*this = other;
			}


			//move assignment operator
			curve& operator=(curve&& other)
			{
				if (this != &other)
				{
					_p = other._p;
					_knot = std::move(other._knot);
					_points = std::move(other._points);
				}
				return *this;
			}
			//move constructor
			curve(curve&& other)
			{
				*this = std::move(other);
			}


			curve(int p, knot<K>& knot, std::vector<T>& points)
			{
				if (p + points.size() + 1 != knot.size()) throw "Illegal parameters: require p + points.size() + 1 = knot.size()";
				_p = p;
				_knot = knot;
				_points = points;
			}

			curve(int p){
				_p = p;
				_knot = geometry::bspline::knot<K>::create_uniform_open(_p, _p + 1);
			}


			T evaluate(K t)
			{
				return geometry::evaluate_curve(_p, _knot.data(), _knot.size(), _points.data(), t);
			}

			void push_back_open(T point)
			{
				_points.push_back(point);
				_knot.push_back_open();
			}
			void push_back_closed(T point)
			{
				_points.push_back(point);
				_knot.push_back_closed();
			}

			void insert_knots(std::vector<K>& k)
			{
				std::vector<K> extra_knots;
				for (int i = 0; i < k.size(); i++)
				{
					if (_knot.minParam() <= k[i] && k[i] <= _knot.maxParam()){
						extra_knots.push_back(k[i]);
					}
				}
				geometry::bspline::knot<K> new_knot(_knot);
				new_knot.insert(extra_knots);

				std::vector<T> new_points(_points.size() + extra_knots.size());

				geometry::olso_insertion(_points.size() - 1, _points.data(), _p + 1, _knot.data(), new_knot.data(), new_knot.size(), new_points.data());
			
				_points = new_points;
				_knot = new_knot;
			}


			void trim_curve(K min_t, K max_t)
			{
				min_t = (min_t > _knot.minParam()) ? min_t : _knot.minParam();
				max_t = (max_t < _knot.maxParam()) ? max_t : _knot.maxParam();

				std::vector<K> extra_knot;
				for (int i = 0; i <= _knot.continuity(min_t); i++){
					extra_knot.push_back(min_t);
				}
				for (int i = 0; i <= _knot.continuity(max_t); i++){
					extra_knot.push_back(max_t);
				}
				this->insert_knots(extra_knot);

				auto k_begin = _knot._data.begin();
				auto p_begin = _points.begin();
				while (*k_begin < min_t) { k_begin++; p_begin++; }
				
				auto k_end = _knot._data.end();
				auto p_end = _points.end();
				while (*(k_end - 1) > max_t) { k_end--; p_end--; }

				_points = std::vector<T>(p_begin, p_end);
				_knot._data = std::vector<K>(k_begin, k_end);

			}



			float length(K minParam, K maxParam)
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
				curve<T,K> c1 = *this;
				c1.trim_curve(min_t, max_t);

				float length1 = 0;
				for (int i = 0; i < c1._points.size() - 1; i++){
					length1 += (c1._points[i] - c1._points[i + 1]).L2norm();
				}
				return length1;
			}


			void set(int i, T point){
				if (i < 0 || i>_n - 1) throw "bspline curve error: require 0 <= i < _n";
				_points[i] = point;
			}
			T get(int i){
				if (i < 0 || i>_n - 1) throw "bspline curve error: require 0 <= i < _n";
				return _points[i];
			}
			T* item(int i){
				if (i < 0 || i>_n - 1) throw "bspline curve error: require 0 <= i < _n";
				return &_points[i];
			}
			inline int nPoints(){
				return _points.size();
			}
			inline int getOrder(){
				return _p;
			}
			//double getKnot(int i){
			//	return knot[i];
			//}
			//double* getKnotVector(){
			//	return knot;
			//}
			//int lKnot(){
			//	return _lknot;
			//}
			inline K minParam(){
				return _knot.minParam();
			}
			inline K maxParam(){
				return _knot.maxParam();
			}







		};



		
		
		//returns true or false if c1 intersects c2.  If there is intersection, parameters of intersect are set.
		template <class T, class K>
		static bool IntersectParam(
			geometry::bspline::curve<vector<T, 2>, K>& c1,
			geometry::bspline::curve<vector<T, 2>, K>& c2,
			T error,
			vector<K, 2>& param,
			bool verbose){

			int max_itr = 100;

			auto d1 = c1;
			auto d2 = c2;

			bool IntersectFound = false;
			bool RefineIntersect = true;
			int i, j;

			for (int itr = 1; itr <= max_itr; itr++){

				if (verbose) { std::cout << "iteration " << itr << std::endl; }

				for (i = 0; i < d1._points.size() - 1; i++){
					for (j = 0; j < d2._points.size() - 1; j++){
						IntersectFound = IntersectionParam(d1._points[i], d1._points[i + 1], d2._points[j], d2._points[j + 1], param, verbose);
						if (IntersectFound)	break;
					}
					if (IntersectFound) break;
				}
				if (IntersectFound){
					double t1 = 0;
					for (int k = i; k < i + d1.getOrder() + 1; k++){
						t1 += param[0] * d1._knot[k] + (1 - param[0]) * d1._knot[k + 1];
					}
					t1 /= (d1.getOrder() + 1);
					double t11[1] = { t1 };
					param[0] = t1;

					double t2 = 0;
					for (int k = j; k < j + d2.getOrder() + 1; k++){
						t2 += param[1] * d2._knot[k] + (1 - param[1]) * d2._knot[k + 1];
					}
					t2 /= (d2.getOrder() + 1);
					double t22[1] = { t2 };
					param[1] = t2;

					auto di1 = d1.evaluate(param[0]);
					auto di2 = d2.evaluate(param[1]);

					d1.insert_knots(std::vector < K > {param[0]});
					d2.insert_knots(std::vector < K > {param[1]});

					if (verbose){
						std::cout << " d1(" << di1[0] << "," << di1[1] << ")";
						std::cout << " d2(" << di2[0] << "," << di2[1] << ")";
						std::cout << std::endl;
						std::cout << " distance = " << sqrt(dot_product(di1, di2)) << std::endl;
					}

					di1 -= di2;
					if (dot_product(di1, di1) < error * error) { return true; }

				}
				else{
					break;
				}
			}
			return false;
		}



	}
}


#endif _GEOMETRY_BSPLINE_CURVE_H_