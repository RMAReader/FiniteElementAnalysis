#ifndef _GEOMETRY_TRIANGLE_H_
#define _GEOMETRY_TRIANGLE_H_

#include "geometry_vector.h"

namespace geometry
{

	template <class T, int N>
	class triangle
	{
	public:
		vector<T, N> p1, p2, p3;

		triangle()
		{
		}

		triangle(vector<T, N>& p1, vector<T, N>& p2, vector<T, N>& p3)
		{
			this->p1 = p1;
			this->p2 = p2;
			this->p3 = p3;
		}

	};


}

#endif _GEOMETRY_TRIANGLE_H_