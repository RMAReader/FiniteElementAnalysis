#ifndef _GEOMETRY_BSPLINE_CURVE_H_
#define _GEOMETRY_BSPLINE_CURVE_H_


#include <vector>
#include "geometry_bspline_knot.h"
#include "bspline_utils.h"

namespace geometry {
	namespace bspline{

		template <class T, class K>
		class curve {
			int _p = 0;
			geometry::bspline::knot<K> _knot;
			std::vector<T> _points;

		public:
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
				if (this != &&other)
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


			curve(int p, std::vector<double>& knot, std::vector<T>& points;)
			{
				if (p + points.size() + 1 != knot.size()) throw "Illegal parameters: require p + points.size() + 1 = knot.size()";
				_p = p;
				_knot = knot;
				_points = points;
			}


			T evaluate(K t)
			{
				return bspline::evaluate_curve(_p, _knot.data(), _knot.size(), _points.data(), t);
			}

			void push_back_open(T point)
			{
				_points.push_back(point);
				_knot.push_back_open();
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

				std::vector<T> new_points(points.size() + extra_knot.size());

				bspline::olso_insertion(points.size() - 1, points.data(), _p + 1, knot.data(), new_knot.data(), new_knot.size(), new_points.data());
			}


			void trim_curve(double min_t, double max_t)
			{
				min_t = maxf(min_t, _knot.minParam());
				max_t = minf(max_t, _knot.maxParam());

				std::vector<K> extra_knot;
				for (int i = 0; i < _knot.continuity(min_t); i++){
					extra_knot.push_back(min_t)
				}
				for (int i = 0; i < _knot.continuity(max_t); i++){
					extra_knot.push_back(max_t)
				}
				this->insert_knots(extra_knot);

				auto k_begin = _knot.begin();
				auto p_begin = _points.begin();
				while (*k_begin < min_t) { k_begin++; p_begin++; }
				
				auto k_end = _knot.end();
				auto p_end = _points.end();
				while (*k_end > max_t) { k_end--; p_end--; }

				_points = std::vector<T>(p_begin, p_end);
				_knot = std::vector<K>(k_begin, k_end);

			}



			T length(K minParam, K maxParam)
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
				curve<T>* c1 = this->trim_curve(min_t, max_t);

				float length1 = 0;
				for (int i = 0; i < c1->nPoints() - 1; i++){
					length1 += distance(c1->get(i), c1->get(i + 1));
				}

				delete c1;
				return length1;
			}


			void set(int i, T point){
				if (i < 0 || i>_n - 1) throw "bspline curve error: require 0 <= i < _n";
				points[i] = point;
			}
			T get(int i){
				if (i < 0 || i>_n - 1) throw "bspline curve error: require 0 <= i < _n";
				return points[i];
			}
			T* item(int i){
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


	}
}


#endif _GEOMETRY_BSPLINE_CURVE_H_