#ifndef _GEOMETRY_BSPLINE_KNOT_H_
#define _GEOMETRY_BSPLINE_KNOT_H_


#include <vector>

namespace geometry
{
	namespace bspline
	{

		template < class T >
		class knot{
			int p;
			std::vector<T> data;

		public:

			knot(){}
			~knot(){}

			//copy assignment
			knot& operator=(const knot& other)
			{
				if (this != &other)
				{
					p = other.p;
					data = other.data;
				}
				return *this;
			}
			//copy constructor
			knot(const knot& other)
			{
				*this = other;
			}


			//move assignment operator
			knot& operator=(knot&& other)
			{
				if (this != &&other)
				{
					p = other.p;
					knot = std::move(other.knot);
				}
				return *this;
			}
			//move constructor
			knot(knot&& other)
			{
				*this = std::move(other);
			}


			static knot& create_uniform_open(int p, int n)
			{
				knot output;
				output.p = p;
				output.data.resize(n);
				for (int i = 0; i < n; i++)
				{
					output[i] = (T)i;
				}
				nromalise();
			}
			static knot& create_uniform_closed(int p, int n)
			{
				knot output;
				output.p = p;
				output.data.resize(n);
				for (int i = 0; i < n; i++){
					output[i] = (T)fmin(fmax((double)(i - p) / (n - p), 0), 1);
				}
				normalise();
			}
			

			inline T minParam()
			{
				return data[p];
			}
			inline T maxParam()
			{
				return data[data.size() - p - 1];
			}
			inline T operator[](int i)
			{
				return data[i];
			}


			void insert(std::vector<T> k)
			{
				data.push_back((data.end(), k.begin(), k.end()));
				std::sort(data.begin(),data.end());
				normalise();
			}
			void push_back_open()
			{
				data = create_uniform_open(p, data.size() + 1);
			}


			//returns continuity of curve at parameter x
			int continuity(K x)
			{
				int c = p;
				for (int i = 0; i < data.size(); i++)
				{
					if (data[i] == x){ c--; }
				}
				return c;
			}
			//returns multiplicity of knot i
			int multiplicity(int i)
			{
				int m = 0;
				for (int j = 0; j < data.size(); j++)
				{
					if (data[i] == data[j]){ m++; }
				}
				return m;
			}

		private:
			
			//rescales knot vector so parameter is in range [0, 1]
			void normalise()
			{
				T alpha = maxParam() - minParam();
				for (int i = 0; i < data.size(); i++)
				{
					data[i] = (data[i] - minParam()) / alpha;
				}
			}

		};

	}
}


#endif _GEOMETRY_BSPLINE_KNOT_H_
