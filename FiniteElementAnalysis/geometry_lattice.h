#ifndef _GEOMETRY_LATTICE_H_
#define _GEOMETRY_LATTICE_H_


#include <vector>

namespace geometry
{


	template <class T>
	class lattice
	{
	public:
		std::vector<T> data;
		int rows, cols;

		lattice(){}

		lattice(lattice& base){
			rows = base.rows;
			cols = base.cols;
			data = base.data;
		}

		lattice(int rows, int cols){
			this->rows = rows;
			this->cols = cols;
			data = std::vector < T >(rows * cols);
			//points.resize(rows * cols);
		}

		lattice(std::vector<T>& p, int rows, int cols){
			if (p.size() == rows * cols){
				this->data = p;
				this->rows = rows;
				this->cols = cols;
			}
		}

		T& GetPoint(int row, int col){
			return data[row * cols + col];
		}
		std::vector<T> GetRow(int row){
			std::vector < T > slice(cols);
			//slice.resize(cols);
			int offset = row * cols;
			for (int i = 0; i < cols; i++){
				slice[i] = data[offset + i];
			}
			return slice;
		}
		std::vector<T> GetCol(int col){
			std::vector < T > slice(rows);
			//slice.resize(cols);
			for (int i = 0; i < rows; i++){
				slice[i] = data[i * cols + col];
			}
			return slice;
		}
		inline void SetPoint(int row, int col, T p){
			data[row * cols + col] = p;
		}
	};


}

#endif _GEOMETRY_LATTICE_H_