#ifndef _GEOMETRY_BSPLINE_KNOT_H_
#define _GEOMETRY_BSPLINE_KNOT_H_


#include <vector>


namespace geometry
{
	namespace bspline
	{

		template < class T >
		class knot{


		public:
			int p;
			std::vector<T> _data;

			knot(){}
			~knot(){}

			//copy assignment
			knot& operator=(const knot& other)
			{
				if (this != &other)
				{
					p = other.p;
					_data = other._data;
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
				if (this != &other)
				{
					p = other.p;
					_data = std::move(other._data);
				}
				return *this;
			}
			//move constructor
			knot(knot&& other)
			{
				*this = std::move(other);
			}


			static knot create_uniform_open(int p, int n)
			{
				knot output;
				output.p = p;
				output._data.resize(n);
				for (int i = 0; i < n; i++)
				{
					output[i] = (T)i;
				}
				output.normalise();
				return output;
			}
			static knot create_uniform_closed(int p, int n)
			{
				knot output = create_uniform_open(p,n);
				output.close_front();
				output.close_back();
				return output;
			}
			

			void close_front()
			{
				for (int i = 0; i < p; i++){
					_data[i] = minParam();
				}
			}
			void close_back()
			{
				for (int i = _data.size() - p; i < _data.size(); i++){
					_data[i] = maxParam();
				}
			}

			inline T minParam()
			{
				return _data[p];
			}
			inline T maxParam()
			{
				return _data[_data.size() - p - 1];
			}
			inline T& operator[](int i)
			{
				return _data[i];
			}


			void insert(std::vector<T> k)
			{
				_data.insert(_data.end(), k.begin(), k.end());
				std::sort(_data.begin(),_data.end());
				normalise();
			}
			void push_back_open()
			{
				*this = create_uniform_open(p, _data.size() + 1);
			}

			void push_back_closed()
			{
				*this = create_uniform_closed(p, _data.size() + 1);
			}

			//returns continuity of curve at parameter x
			int continuity(T x)
			{
				int c = p;
				for (int i = 0; i < _data.size(); i++)
				{
					if (_data[i] == x){ c--; }
				}
				return c;
			}
			//returns multiplicity of knot i
			int multiplicity(int i)
			{
				int m = 0;
				for (int j = 0; j < _data.size(); j++)
				{
					if (_data[i] == _data[j]){ m++; }
				}
				return m;
			}

			inline int size()
			{
				return _data.size();
			}
			inline T* data()
			{
				return _data.data();
			}


		private:
			
			//rescales knot vector so parameter is in range [0, 1]
			void normalise()
			{
				T alpha = maxParam() - minParam();
				T beta = minParam();
				for (int i = 0; i < _data.size(); i++)
				{
					_data[i] = (_data[i] - beta) / alpha;
				}
			}

		};

	}
}


#endif _GEOMETRY_BSPLINE_KNOT_H_
