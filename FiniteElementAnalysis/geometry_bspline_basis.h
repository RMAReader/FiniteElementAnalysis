#ifndef _GEOMETRY_BSPLINE_BASIS_H_
#define _GEOMETRY_BSPLINE_BASIS_H_


namespace geometry {
	namespace bspline{

		template < class T_k, class T_x >
		class basis
		{
			T_k _knot;
			int _p;

		public:
			basis()
			{
				
			}

			//copy assignment
			basis& operator=(const basis& other)
			{
				if (this != &other)
				{
					_knot = other._knot;
					_p = other._p
				}
				return *this;
			}
			//copy constructor
			basis(const basis& other)
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


			public T_x evaluate(int i, int p, double x){//, double *knot, int lknot, double x){

				if (i < 0) throw "deboor error: require i >= 0";
				
				if (i + p + 1 >= knot.size()) throw "deboor error: require i + p + 1 < knot.size()";

				if (x == knot[knot.size() - p - 1]){ x -= x * DBL_EPSILON; }

				if (p == 0){
					if (knot[i] <= x && x < knot[i + 1]){
						return 1;
					}
					else{
						return 0;
					}
				}
				double d1 = knot[i + p] - knot[i];
				double d2 = knot[i + p + 1] - knot[i + 1];

				double val1 = 0;
				double val2 = 0;
				if (d1 > 0){
					val1 = (x - knot[i]) / d1 * deboor_value(i, p - 1, knot, lknot, x);
				}
				if (d2 > 0) {
					val2 = (knot[i + p + 1] - x) / d2 * deboor_value(i + 1, p - 1, knot, lknot, x);
				}
				return val1 + val2;
			}



		};



	}
}


#endif _GEOMETRY_BSPLINE_BASIS_H_