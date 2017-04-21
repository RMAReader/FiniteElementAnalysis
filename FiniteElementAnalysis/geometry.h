#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

#include <iostream>

#include "geometry_vector.h"
#include "geometry_matrix.h"
#include "geometry_bspline_basis.h"
#include "geometry_bspline_curve.h"
#include "geometry_bspline_knot.h"
#include "geometry_bspline_algorithms.h"
#include "geometry_mesh.h"
#include "geometry_bspline_surface.h"
#include "geometry_circle.h"

#define geoKNOT geometry::bspline::knot<double>
#define geoCURVE2F geometry::bspline::curve<geometry::vector<float, 2>, double>
#define geoCURVE3F geometry::bspline::curve<geometry::vector<float, 3>, double>
#define geoVEC2F geometry::vector<float, 2>
#define geoVEC3F geometry::vector<float, 3>
#define geoVEC2D geometry::vector<double, 2>
#define geoVEC3D geometry::vector<double, 3>
#define geoSURFACE3F geometry::bspline::surface<geometry::vector<float, 3>, double, double>
#define geoLATTICE3F geometry::lattice<geometry::vector<float,3>>
#define geoCIRCLE2F geometry::circle<float,2>

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

	//returns true if point lies on or within circle
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
	static bool interior_triangle(vector<T, 2>& c, T r, vector<T, 2>& p1, vector<T, 2>& p2, vector<T, 2>& p3){

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

	template <class T>
	static T min_height_sphere_on_point(T x, T y, T r, vector<T, 3>& p)
	{
		vector<T, 2> c2; c2[0] = x; c2[1] = y;
		vector<T, 2> p2; p2 << p;

		T distance = (p2 - c2).L2norm();
		
		if (distance <= r)
		{
			return p[2] - r + sqrt(r * r - distance * distance);
		}
		else
		{
			return DBL_MIN;
		}
	}



	template <class T>
	static bool min_height_sphere_on_triangle(T x, T y, T r, vector<T, 3>& p1, vector<T, 3>& p2, vector<T, 3>& p3, T& h)
	{
		
		//1. tool tip touches triangle surface
		vector<T, 3> q1;
		vector<T, 3> q2 = p2 - p1;
		vector<T, 3> q3 = p3 - p1;

		vector<T, 3> n = geometry::cross_product(q2, q3);

		T det = n[2];
		if (n[2] < 0) n *= -1;

		n.normalise();

		q1[0] = x - n[0] - p1[0];
		q1[1] = y - n[1] - p1[1];

		T alpha = (q3[1] * q1[0] - q3[0] * q1[1]) / det;
		T beta = (-q2[1] * q1[0] + q2[0] * q1[1]) / det;

		if (0 <= alpha && alpha <= 1 && 0 <= beta && beta <= 1 && alpha + beta <= 1)
		{
			h= p1[2] + (r - n[0] * (x - p1[0]) - n[1] * (y - p1[1])) / n[2];
			return true;
		}

		//2. tool tip touches triangle edge
		T z;
		bool contact = false;

		geometry::mesh_edge<float> e;
		e.initialise(p1, p2);
		if (e.min_height_sphere(x, y, r, z))
		{
			if (contact){ h = fmaxf(h, z); }
			else{ h = z; }
			contact = true;
		}
		e.initialise(p2, p3);
		if (e.min_height_sphere(x, y, r, z))
		{
			if (contact){h = fmaxf(h, z);}
			else{ h = z; }
			contact = true;
		}
		e.initialise(p3, p1);
		if (e.min_height_sphere(x, y, r, z))
		{
			if (contact){ h = fmaxf(h, z); }
			else{ h = z; }
			contact = true;
		}
		if (contact)return true;

		//3. tool tip touches triangle vertex
		T  x2, y2,d,r2;
		r2 = r*r;

		x2 = x - p1[0]; x2 *= x2;
		y2 = y - p1[1]; y2 *= y2;
		d = r2 - x2 - y2;
		if (d >= 0)
		{
			h= p1[2] + sqrtf(d);
			return true;
		}
		x2 = x - p2[0]; x2 *= x2;
		y2 = y - p2[1]; y2 *= y2;
		d = r2 - x2 - y2;
		if (d >= 0)
		{
			h= p2[2] + sqrtf(d);
			return true;
		}
		x2 = x - p3[0]; x2 *= x2;
		y2 = y - p3[1]; y2 *= y2;
		d = r2 - x2 - y2;
		if (d >= 0)
		{
			h= p3[2] + sqrtf(d);
			return true;
		}
		return false;
	}

	template <class T>
	static std::vector<vector<T, 3>> intersection_triangle_plane_xz(T y, vector<T, 3>& p1, vector<T, 3>& p2, vector<T, 3>& p3)
	{
		std::vector<vector<T, 3>> segment(2);
		vector<T, 3> intersect;

		if (intersection_line_plane_xz(y, p1, p2, intersect))
		{
			segment.push_back(intersect);
		}
		if (intersection_line_plane_xz(y, p2, p3, intersect))
		{
			segment.push_back(intersect);
		}
		if (intersection_line_plane_xz(y, p3, p1, intersect))
		{
			segment.push_back(intersect);
		}
		return segment;
	}

	template <class T>
	static bool intersection_line_plane_xz(T y, vector<T, 3>& p1, vector<T, 3>& p2, vector<T, 3>& intersect)
	{
		if (p1[1] <= y && p2[1] > y)
		{
			T s = (y - p1[1]) / (p2[1] - p1[1]);
			intersect = s * p1 + (1 - s) * p2;
			return true;
		}
		else if (p2[1] <= y && p1[1] > y)
		{
			T s = (p1[1] - y) / (p1[1] - p2[1]);
			intersect = s * p1 + (1 - s) * p2;
			return true;
		}
		return false;
	}


}


#endif _GEOMETRY_H_
