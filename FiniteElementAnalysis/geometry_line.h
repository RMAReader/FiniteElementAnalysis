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

		inline void reverse()
		{
			std::swap(p1, p2);
		}
	};


	template <class T, int N, int i>
	inline void sort_line(line<T, N>& L1)
	{
		if (L1.p1[i] > L1.p2[i]){ L1.reverse(); }
	}
	

	template <class T, int N, int i>
	class sort_ascending_p1
	{
	public:
		inline bool operator() (const line<T, N>& L1, const line<T, N>& L2)
		{
			return (L1.p1[i] < L2.p1[i]);
		}
	};

	template <class T, int N>
	class sort_ascending_p2
	{
	private:
		int i;
	public:
		sort_ascending_p2(int i){ this->i = i; }
		inline bool operator() (const line<T, N>& L1, const line<T, N>& L2)
		{
			return (L1.p2[i] < L2.p2[i]);
		}
	};

}

#endif _GEOMETRY_LINE_H_