#ifndef _GEOMETRY_RECTANGLE_H_
#define _GEOMETRY_RECTANGLE_H_

#include "geometry_vector.h"

namespace geometry
{

	template <class T, int N>
	class rectangle
	{
	public:
		vector<T, N> minimum;
		vector<T, N> maximum;

		rectangle()
		{
		}

		rectangle(vector<T, N> minimum, vector<T, N> maximum)
		{
			this->minimum = minimum;
			this->maximum = maximum;
		}

	};


}

#endif _GEOMETRY_RECTANGLE_H_