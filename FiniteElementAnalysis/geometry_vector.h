#ifndef _GEOMETRY_VECTOR_H_
#define _GEOMETRY_VECTOR_H_

#include <array>


namespace geometry
{

	template <class T, int _N>
	class vector
	{
	private:
		
		std::array<T, _N> data;

	public:

		vector()
		{
			for (int i = 0; i < _N; i++){
				this->data[i] = 0;
			}
		}

		//copy assignment
		vector& vector::operator=(const vector& other)
		{
			if (this != &other)
			{
				for (int i = 0; i < _N; i++){
					this->data[i] = other.data[i];
				}
			}
			return *this;
		}
		//copy constructor
		vector(const vector& other)
		{
			*this = other;
		}

		vector(const std::array<T, _N>& data)
		{
			this->data = data;
		}



		////move assignment operator
		//vector& operator=(vector&& other)
		//{
		//	if (this != &&other)
		//	{
		//		data = std::move(other.data);
		//	}
		//	return *this;
		//}
		////move constructor
		//vector::vector(vector&& other)
		//{
		//	*this = std::move(other);
		//}


		//vector::vector(T* data, int s)
		//{
		//	for (int i = 0; i < fminf(s,N); i++){
		//		this->data[i] = other.data[i];
		//	}
		//	for (int i = s; i < N; i++){
		//		this->data[i] = 0;
		//	}
		//}


//		virtual ~vector(){}

		//copy assignment
		void operator=(T v)
		{
			for (int i = 0; i < _N; i++){
				this->data[i] = v;
			}
		}
		void operator +=(vector& v)
		{
			for (int i = 0; i < _N; i++){
				this->data[i] += v.data[i];
			}
		}
		void operator +=(T v)
		{
			for (int i = 0; i < _N; i++){
				this->data[i] += v;
			}
		}
		void operator -=(vector& v)
		{
			for (int i = 0; i < _N; i++){
				this->data[i] -= v.data[i];
			}
		}
		void operator -=(T v)
		{
			for (int i = 0; i < _N; i++){
				this->data[i] -= v;
			}
		}
		void operator *=(vector& v)
		{
			for (int i = 0; i < _N; i++){
				this->data[i] *= v.data[i];
			}
		}
		void operator *=(T v)
		{
			for (int i = 0; i < _N; i++){
				this->data[i] *= v;
			}
		}
		void operator /=(vector& v)
		{
			for (int i = 0; i < _N; i++){
				this->data[i] /= v.data[i];
			}
		}
		void operator /=(T v)
		{
			for (int i = 0; i < _N; i++){
				this->data[i] /= v;
			}
		}

		inline T& operator [](int i)
		{
			if (0 > i || i >= _N) throw "Argument out of range";
			return data[i];
		}

		inline T L1norm()
		{
			T out = 0;
			for (int i = 0; i < _N; i++)
			{
				out += abs(data[i]);
			}
			return out;
		}
		
		inline T L2norm()
		{
			T out = 0;
			for (int i = 0; i < _N; i++)
			{
				out += data[i] * data[i];
			}
			return sqrt(out);
		}

		void normalise()
		{
			T s = 1 / this->L2norm();
			for (int i = 0; i < _N; i++)
			{
				data[i] *= s;
			}
		}

		inline int size()
		{
			return _N;
		}
	};



	template <class T, int N, int M>
	static void operator<<(vector<T, N>& left, vector<T, M>& right)
	{
		int size = (N < M) ? N : M;
		for (int i = 0; i < size; i++)
		{
			left[i] = right[i];
		}
	}


	template <class T, int N>
	static vector<T, N> operator+(vector<T, N>& a, vector<T, N>& b)
	{
		vector<T, N> out = a;
		out += b;
		return out;
	}
	template <class T, int N>
	static vector<T, N> operator+(T a, vector<T, N>& b)
	{
		vector<T, N> out = b;
		out += a;
		return out;
	}
	template <class T, int N>
	static vector<T, N> operator+(vector<T, N>& a, T b)
	{
		vector<T, N> out = a;
		out += b;
		return out;
	}

	template <class T, int N>
	static vector<T, N> operator-(vector<T, N>& a, vector<T, N>& b)
	{
		vector<T, N> out = a;
		out -= b;
		return out;
	}

	template <class T, int N>
	static vector<T, N> operator*(vector<T, N>& a, vector<T, N>& b)
	{
		vector<T, N> out = a;
		out *= b;
		return out;
	}
	template <class T, int N>
	static vector<T, N> operator*(vector<T, N>& a, T b)
	{
		vector<T, N> out = a;
		out *= b;
		return out;
	}
	template <class T, int N>
	static vector<T, N> operator*(T a, vector<T, N>& b)
	{
		vector<T, N> out = b;
		out *= a;
		return out;
	}

	template <class T, int N>
	static vector<T, N> operator/(vector<T, N>& a, vector<T, N>& b)
	{
		vector<T, N> out = a;
		out /= b;
		return out;
	}



}

#endif _GEOMETRY_VECTOR_H_