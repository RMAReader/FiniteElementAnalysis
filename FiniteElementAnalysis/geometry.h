#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

#include "geometry_vector.h"
#include "geometry_matrix.h"
#include "geometry_bspline_basis.h"

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





}


#endif _GEOMETRY_H_
