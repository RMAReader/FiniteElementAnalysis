#ifndef _GEOMETRY_LINE_H_
#define _GEOMETRY_LINE_H_

#include "geometry_vector.h"

namespace geometry
{

	template <class T, int N>
	class line
	{
	public:
		vector<T, N> p1, p2;

		line()
		{
		}

		line(vector<T, N>& p1, vector<T, N>& p2)
		{
			this->p1 = p1;
			this->p2 = p2;
		}

	};


}

#endif _GEOMETRY_LINE_H_