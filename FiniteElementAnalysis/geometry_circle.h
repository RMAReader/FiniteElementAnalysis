#ifndef _GEOMETRY_CIRCLE_H_
#define _GEOMETRY_CIRCLE_H_

#include "geometry_vector.h"

namespace geometry
{

	template <class T, int N>
	class circle
	{
	public:
		T radius;
		vector<T, N> centre;

		circle(vector<T, N> centre, T radius)
		{
			this->centre = centre;
			this->radius = radius;
		}

	};


}

#endif _GEOMETRY_CIRCLE_H_