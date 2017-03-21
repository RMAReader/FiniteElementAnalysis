#ifndef _GEOMETRY_MATRIX_H_
#define _GEOMETRY_MATRIX_H_

#include "geometry_vector.h"

namespace geometry
{

	template <class T, int _ROWS, int _COLS>
	class matrix
	{
	private:
		//int N = _N;
		T data[_ROWS*_COLS];


	public:

		matrix()
		{
			for (int i = 0; i < _ROWS*_COLS; i++){
				this->data[i]= 0;
			}
		}

		//copy assignment
		matrix& matrix::operator=(const matrix& other)
		{
			if (this != &other)
			{
				for (int i = 0; i < _ROWS*_COLS; i++){
					this->data[i] = other.data[i];
				}
			}
			return *this;
		}
		//copy constructor
		matrix(const matrix& other)
		{
			*this = other;
		}


		////move assignment operator
		//matrix& matrix::operator=(matrix&& other)
		//{
		//	if (this != &other)
		//	{
		//		for (int i = 0; i < N; i++){
		//			this->data[i] = other.data[i];
		//		}
		//	}
		//	return *this;
		//}
		////move constructor
		//matrix::matrix(matrix&& other)
		//{
		//	*this = std::move(other);
		//}


		//matrix::matrix(T* data, int s)
		//{
		//	for (int i = 0; i < fminf(s,N); i++){
		//		this->data[i] = other.data[i];
		//	}
		//	for (int i = s; i < N; i++){
		//		this->data[i] = 0;
		//	}
		//}


		virtual ~matrix(){}


		void operator +=(matrix& v)
		{
			for (int i = 0; i < _ROWS*_COLS; i++){
				this->data[i] += v.data[i];
			}
		}
		void operator -=(matrix& v)
		{
			for (int i = 0; i < _ROWS*_COLS; i++){
				this->data[i] -= v.data[i];
			}
		}
		void operator *=(matrix& v)
		{
			for (int i = 0; i <_ROWS*_COLS; i++){
				this->data[i] *= v.data[i];
			}
		}
		void operator /=(matrix& v)
		{
			for (int i = 0; i < _ROWS*_COLS; i++){
				this->data[i] /= v.data[i];
			}
		}

		inline T& operator ()(int i,int j)
		{
			if (0 > i || i >= _ROWS) throw "Argument out of range";
			if (0 > j || j >= _COLS) throw "Argument out of range";
			return data[i * _COLS + j];
		}

		inline int size()
		{
			return _ROWS*_COLS;
		}
		inline int rows()
		{
			return _ROWS;
		}
		inline int cols()
		{
			return _COLS;
		}



	};


	
	template <class T, int _ROWS,int  _COLS>
	static matrix<T, _ROWS, _COLS> operator+(matrix<T, _ROWS, _COLS>& a, matrix<T, _ROWS, _COLS>& b)
	{
		vector<T, N> out = a;
		out += b;
		return out;
	}

	template <class T, int _ROWS, int _COLS>
	static matrix<T, _ROWS, _COLS> operator-(matrix<T, _ROWS, _COLS>& a, matrix<T, _ROWS, _COLS>& b)
	{
		vector<T, N> out = a;
		out -= b;
		return out;
	}

	template <class T,int  _ROWS, int _COLS>
	static matrix<T, _ROWS, _COLS> operator*(matrix<T, _ROWS, _COLS>& a, matrix<T, _ROWS, _COLS>& b)
	{
		vector<T, N> out = a;
		out *= b;
		return out;
	}

	template <class T, int _ROWS, int _COLS>
	static matrix<T, _ROWS, _COLS> operator/(matrix<T, _ROWS, _COLS>& a, matrix<T, _ROWS, _COLS>& b)
	{
		vector<T, N> out = a;
		out /= b;
		return out;
	}


	
	

}

#endif _GEOMETRY_MATRIX_H_
