#ifndef _GEOMETRY_VECTOR_H_
#define _GEOMETRY_VECTOR_H_


namespace geometry
{

	template <class T, int _N>
	class vector
	{
	private:
		T data[_N];
		

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


		////move assignment operator
		//vector& vector::operator=(vector&& other)
		//{
		//	if (this != &&other)
		//	{
		//		for (int i = 0; i < N; i++){
		//			this->data[i] = other.data[i];
		//		}
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


		virtual ~vector(){}


		void operator +=(vector& v)
		{
			for (int i = 0; i < _N; i++){
				this->data[i] += v.data[i];
			}
		}
		void operator -=(vector& v)
		{
			for (int i = 0; i < _N; i++){
				this->data[i] -= v.data[i];
			}
		}
		void operator *=(vector& v)
		{
			for (int i = 0; i < _N; i++){
				this->data[i] *= v.data[i];
			}
		}
		void operator /=(vector& v)
		{
			for (int i = 0; i < _N; i++){
				this->data[i] /= v.data[i];
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



	template <class T, int N>
	static vector<T, N> operator+(vector<T, N>& a, vector<T, N>& b)
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
	static vector<T, N> operator/(vector<T, N>& a, vector<T, N>& b)
	{
		vector<T, N> out = a;
		out /= b;
		return out;
	}

}

#endif _GEOMETRY_VECTOR_H_