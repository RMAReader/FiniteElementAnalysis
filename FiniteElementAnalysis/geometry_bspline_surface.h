#ifndef _GEOMETRY_BSPLINE_SURFACE_H_
#define _GEOMETRY_BSPLINE_SURFACE_H_


#include <vector>
//#include <algorithm>  

#include "geometry_bspline_algorithms.h"
#include "geometry_bspline_curve.h"
#include "geometry_bspline_knot.h"
#include "geometry_lattice.h"


namespace geometry {
	namespace bspline{


		template <class T, class U, class V>
		class surface
		{
		public:
			int _p;
			int _q;
			lattice<T> _points;
			knot<U> _knotx;
			knot<V> _knoty;

			surface(){}

			//copy assignment
			surface& operator=(const surface& other)
			{
				if (this != &other)
				{
					_p = other._p;
					_q = other._q;
					_knotx = other._knotx;
					_knoty = other._knoty;
					_points = other._points;
				}
				return *this;
			}
			//copy constructor
			surface(const surface& other)
			{
				*this = other;
			}


			//move assignment operator
			surface& operator=(surface&& other)
			{
				if (this != &other)
				{
					_p = other._p;
					_q = other._q;
					_knotx = std::move(other._knotx);
					_knoty = std::move(other._knoty);
					_points = std::move(other._points);
				}
				return *this;
			}
			//move constructor
			surface(surface&& other)
			{
				*this = std::move(other);
			}



			surface(int p, int q, knot<U>& knotx, knot<V>& knoty, lattice<T>& points)
			{
				this->_p = p;
				this->_q = q;
				this->_points = points;
				this->_knotx = knotx;
				this->_knoty = knoty;
			}

			curve<T,U> curve_x(double v){
				std::vector<T> c(points.rows);
				for (int j = 0; j < points.rows; j++){
					c[j] = evaluate_curve(_q, _knoty.data(), _knoty.size(), _points.GetCRow(j).data(), v);
				}
				return curve<T,U>(_p, _knotx, c);
			}

			curve<T,V> curve_y(double u){
				std::vector<T> c(_points.cols);
				for (int j = 0; j < _points.cols; j++){
					c[j] = evaluate_curve(_p, _knotx.data(), _knotx.size(), _points.GetCol(j).data(), u);
				}
				return curve<T,V>(_q, _knoty, c);
			}

			T evaluate(double u, double v){
				std::vector<T> c(points.cols);
				for (int j = 0; j < points.cols; j++){
					c[j] = evaluate_curve(p, *(knotx.data()), knotx.size(), *(points.GetCol(j).data()), u);
				}
				return evaluate_curve(q, *(knoty.data()), knoty.size(), *(c.data()), v);
			}




		};


	}
}


#endif _GEOMETRY_BSPLINE_SURFACE_H_