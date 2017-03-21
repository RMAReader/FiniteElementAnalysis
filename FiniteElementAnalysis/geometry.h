#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

#include "geometry_vector.h"
#include "geometry_matrix.h"
#include "geometry_bspline_basis.h"
#include "geometry_bspline_curve.h"
#include "geometry_bspline_knot.h"
#include "geometry_bspline_algorithms.h"
#include "geometry_mesh.h"
#include "geometry_bspline_surface.h"

#define geoCURVE2F geometry::bspline::curve<geometry::vector<float, 2>, double>
#define geoCURVE3F geometry::bspline::curve<geometry::vector<float, 3>, double>
#define geoVEC2F geometry::vector<float, 2>
#define geoVEC3F geometry::vector<float, 3>
#define geoVEC2D geometry::vector<double, 2>
#define geoVEC3D geometry::vector<double, 3>
#define geoSURFACE3F geometry::bspline::surface<geometry::vector<float, 3>, double, double>
#define geoLATTICE3F geometry::lattice<geometry::vector<float,3>>

namespace geometry
{

	template <class T, int ROWS, int COLS>
	static vector<T, ROWS>& product(matrix<T, ROWS, COLS>& m, vector<T, COLS>& v)
	{
		vector<T, ROWS> out;
		for (int i = 0; i < ROWS; i++){
			for (int j = 0; j < COLS; j++){
				out[i] += m(i, j) * v[j];
			}
		}
		return out;
	}
	
	template <class T, int ROWS, int COLS>
	static vector<T, COLS>& product(vector<T, ROWS>& v, matrix<T, ROWS, COLS>& m)
	{
		vector<T, COLS> out;
		for (int j = 0; j < COLS; j++){
			for (int i = 0; i < ROWS; i++){
				out[j] += v[i] * m(i, j);
			}
		}
		return out;
	}

	template <class T, int N>
	static T dot_product(vector<T, N>& a, vector<T, N>& b)
	{
		T out = 0;
		for (int i = 0; i < N; i++)
		{
			out += a[i] * b[i];
		}
		return out;
	}


	template <class T>
	static vector<T, 3> cross_product(vector<T, 3>& a, vector<T, 3>& b)
	{
		vector<T, 3> out;
		out[0] = a[1] * b[2] - a[2] * b[1];
		out[1] = a[2] * b[0] - a[0] * b[2];
		out[2] = a[0] * b[1] - a[1] * b[0];
		return out;
	}

	template <class T, class K>
	static bool IntersectionParam(vector<T,2>& p1, vector<T,2>& p2, vector<T,2>& q1, vector<T,2>& q2, vector<K,2>& param, bool verbose){
		float a = p1[0] - p2[0];
		float b = p1[1] - p2[1];
		float c = q1[0] - q2[0];
		float d = q1[1] - q2[1];

		double det = a * d - b * c;
		if (det == 0){ return false; }

		float e = q2[0] - p2[0];
		float f = q2[1] - p2[1];

		double alpha = (d * e - c * f) / det;
		double beta = (b * e - a * f) / det;

		param[0] = (K)alpha;
		param[1] = (K)beta;

		if (alpha >= 0 && alpha <= 1 && beta >= 0 && beta <= 1){

			if (verbose){
				std::cout << "intersection: " << std::endl;
				std::cout << "  p1=(" << p1[0] << "," << p1[1] << ") p2=(" << p2[0] << "," << p2[1] << ")" << std::endl;
				std::cout << "  q1=(" << q1[0] << "," << q1[1] << ") q2=(" << q2[0] << "," << q2[1] << ")" << std::endl;
				std::cout << "  alpha=" << alpha << ", beta=" << beta << std::endl;
			}
			return true;
		}
		else{
			return false;
		}
	}

	template <class T>
	static bool interior_circle_point(vector<T, 2>& c, T r, vector<T, 2> p1){
		p1 -= c;
		if (geometry::dot_product(p1,p1) < r * r){ return true; }
		return false;
	}


	//returns true if circle overlaps p1->p2 edge of triangle such that could overlap interior of triangle [p1, p2, p3]
	template <class T>
	static bool interior_triangle_edge(vector<T, 2> c, T r, vector<T, 2> p1, vector<T, 2> p2, vector<T, 2> p3){

		//1. translate points so p1 is origin
		p2 -= p1;
		p3 -= p1;
		c -= p1;
		p1 -= p1;

		//2. rotate points about origin so p2 lies on positive x axis
		vector<T, 2> q2;
		q2[0] = p2.L2norm();
		q2[1] = 0;

		vector<T, 2> d;
		d[0] = (c[0] * p3[0] + c[1] * p3[1]) / q2[0]; 
		d[1] = (c[0] * p3[1] - c[1] * p3[0]) / q2[0];
		

		if (d[0] + r < 0){ return false; }
		if (d[0] - r > q2[0]){ return false; }

		vector<T, 2> q3;
		q3[0] = (p2[0] * p3[0] + p2[1] * p3[1]) / q2[0];
		q3[1] = (p2[0] * p3[1] - p2[1] * p3[0]) / q2[0];

		if (d[1] + r < fminf(0, q3[1])){ return false; }
		if (d[1] - r > fmaxf(0, q3[1])){ return false; }

		return true;
	}

	//returns true if circle and triangle overlap
	template <class T>
	static bool interior_triangle(vector<T, 2>& c, float r, vector<T, 2>& p1, vector<T, 2>& p2, vector<T, 2>& p3){

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


}


#endif _GEOMETRY_H_
