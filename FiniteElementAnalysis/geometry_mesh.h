#ifndef _GEOMETRY_MESH_H_
#define _GEOMETRY_MESH_H_

#include <vector>
#include <unordered_set>
#include "geometry_vector.h"
#include "geometry_lattice.h"
#include "geometry_bspline_surface.h"
#include "geometry_line.h"
#include "geometry_triangle.h"


namespace geometry
{


	class mesh3f
	{
	public:

		std::vector<geometry::vector<float,3>> points;
		std::vector<std::vector<int>> elements;

		lattice<std::vector<int>> nodes;
		std::vector<float> node_x;
		std::vector<float> node_y;


		mesh3f();
		mesh3f(geometry::bspline::surface<geometry::vector<float, 3 >,double, double>& surface, int nx, int ny);


		void add_triangle(triangle<float,3>&);

		//inline triangle3f get_triangle(int index);
		geometry::vector<float, 3> get_vertex(int element, int v);
		//std::vector<triangle3f>* triangles();

		void get_mesh_range(geometry::vector<float, 3>& min, geometry::vector<float, 3>& max);
		inline void get_element_range(int i, geometry::vector<float, 3>& min, geometry::vector<float, 3>& max);
		void build_nodes(int nx, int ny);
		void build_grid_z(float, float);
	};


	template <class T>
	class mesh_triangle
	{
	private:

	public:
		vector<T, 3> p1;
		vector<T, 3> p2;
		vector<T, 3> p3;

		mesh_triangle(vector<T, 3>& v1, vector<T, 3>& v2, vector<T, 3>& v3)
		{
			p1 = v1;
			p2 = v2;
			p3 = v3;
		}
			
		bool min_height_sphere(T x, T y, T r, T& h)
		{
			//1. tool tip touches triangle surface
			vector<T, 3> q1;
			vector<T, 3> q2 = p2 - p1;
			vector<T, 3> q3 = p3 - p1;

			vector<T, 3> n = geometry::cross_product(q2, q3);

			T det = n[2];
			if (n[2] < 0) n *= -1;

			n.normalise();

			q1[0] = x - n[0] - p1[0];
			q1[1] = y - n[1] - p1[1];

			T alpha = (q3[1] * q1[0] - q3[0] * q1[1]) / det;
			T beta = (-q2[1] * q1[0] + q2[0] * q1[1]) / det;

			if (0 <= alpha && alpha <= 1 && 0 <= beta && beta <= 1 && alpha + beta <= 1)
			{
				h = p1[2] + alpha * q2[2] + beta * q3[2] + r * (n[2] - 1);
				return true;
			}
		}
	};

	template <class T>
	class mesh_edge
	{
	private:
		T l, m, m0, m1, n, ux, uy, uz, a;
		vector<T, 2> q;
		vector<T, 3> t;
	public:
		vector<T, 3> p1;
		vector<T, 3> p2;

		mesh_edge(vector<T, 3>& v1, vector<T, 3>& v2)
		{
			p1 = v1;
			p2 = v2;

			t = p2 - p1;
			l = t[0] * t[0] + t[1] * t[1] + t[2] * t[2];
			m0 = -t[0] / l;
			m1 = -t[1] / l;
			m = -t[2] / l;
			ux = m * t[0];
			uy = m * t[1];
			uz = m * t[2] + 1;
			a = ux*ux + uy*uy + uz*uz;
			ux /= a;
			uy /= a;
			uz /= a;
			a *= 2;
		}

		bool min_height_sphere(T x, T y, T r, T& h)
		{
			T q0, q1, vx, vy, vz, b, c, d, r2, z, alpha;

			q0 = p1[0] - x;
			q1 = p1[1] - y;

			n = q0 * m0 + q1 * m1;

			vx = n * t[0] + q0;
			vy = n * t[1] + q1;
			vz = n * t[2];

			b = ux*vx + uy*vy + uz*vz;
			c = vx*vx + vy*vy + vz*vz - r*r;
			
			d = b*b - 2 * a*c;
			bool contact = false;
			if (d >= 0)
			{
				d = sqrtf(d) / a;
				z = -b + d;
				alpha = m * z + n;

				if (0 <= alpha && alpha <= 1){
					z = p1[2] - z;
					h = z - r;
					contact = true;
				}
				
				z = -b - d;
				alpha = m * z + n;
				if (0 <= alpha && alpha <= 1){
					z = p1[2] - z;
					if (contact)
					{
						h = fmaxf(h, z - r);
					}
					else{
						h = z - r;
						contact = true;
					}
				}
			}
			return contact;
		}

	};



	class mesh3f_region
	{
	private:
		mesh3f* mesh;
		
		int min_x_index;
		int min_y_index;
		int max_x_index;
		int max_y_index;

		float min_x;
		float min_y;
		float max_x;
		float max_y;
		int nx;
		int ny;

		std::vector<mesh_triangle<float>> all_triangles;
		std::vector<mesh_edge<float>> all_edges;
		std::vector<geometry::vector<float, 3>> all_points;


		lattice<std::vector< int >> nodes_points;
		lattice<std::vector< int >> nodes_edges;
		lattice<std::vector< int >> nodes_triangles;


	public:

		mesh3f_region();
		mesh3f_region(mesh3f* mesh, int nx, int ny);

		void build_nodes(mesh3f* mesh, int nx, int ny);
		void set_region(float min_x, float max_x, float min_y, float max_y);

		std::vector<int> distinct_triangles;
		std::vector<int> distinct_edges;
		std::vector<int> distinct_points;

		typedef std::vector<int>::iterator triangle_iterator;
		triangle_iterator begin_triangle() { return distinct_triangles.begin(); }
		triangle_iterator end_triangle()   { return distinct_triangles.end(); }

		typedef std::vector<int>::iterator edge_iterator;
		edge_iterator begin_edge() { return distinct_edges.begin(); }
		edge_iterator end_edge()   { return distinct_edges.end(); }

		typedef std::vector<int >::iterator point_iterator;
		point_iterator begin_point() { return distinct_points.begin(); }
		point_iterator end_point()   { return distinct_points.end(); }

		inline mesh_triangle<float> get_triangle(int i){ return all_triangles[i]; }
		inline mesh_edge<float> get_edge(int i){ return all_edges[i]; }
		inline geometry::vector<float, 3> get_point(int i){ return mesh->points[i]; }
	};





	////template <class T>
	//class triangle_iterator : public std::iterator<std::random_access_iterator_tag,
	//	geometry::triangle<float, 3>,
	//	ptrdiff_t,
	//	geometry::triangle<float, 3>*,
	//	geometry::triangle<float, 3>&>
	//{
	//public:

	//	triangle_iterator(mesh3f_region* region = nullptr){ this->region = region; }
	//	triangle_iterator(const triangle_iterator& rawIterator) = default;
	//	~triangle_iterator(){}

	//	triangle_iterator& operator=(const triangle_iterator& rawIterator) = default;
	//	triangle_iterator& operator=(geometry::triangle<float, 3>* ptr){ m_ptr = ptr; return (*this); }

	//	operator bool()const
	//	{
	//		if (m_ptr)
	//			return true;
	//		else
	//			return false;
	//	}

	//	bool operator==(const triangle_iterator& rawIterator)const{ return (m_ptr == rawIterator.getConstPtr()); }
	//	bool operator!=(const triangle_iterator& rawIterator)const{ return (m_ptr != rawIterator.getConstPtr()); }

	//	triangle_iterator& operator+=(const ptrdiff_t& movement){ m_ptr += movement; return (*this); }
	//	triangle_iterator& operator-=(const ptrdiff_t& movement){ m_ptr -= movement; return (*this); }
	//	triangle_iterator& operator++(){ ++m_ptr; return (*this); }
	//	triangle_iterator& operator--(){ --m_ptr; return (*this); }
	//	//triangle_iterator  operator++(ptrdiff_t){ auto temp(*this); ++m_ptr; return temp; }
	//	//triangle_iterator  operator--(ptrdiff_t){ auto temp(*this); --m_ptr; return temp; }
	//	triangle_iterator  operator+(const ptrdiff_t& movement){ auto oldPtr = m_ptr; m_ptr += movement; auto temp(*this); m_ptr = oldPtr; return temp; }
	//	triangle_iterator  operator-(const ptrdiff_t& movement){ auto oldPtr = m_ptr; m_ptr -= movement; auto temp(*this); m_ptr = oldPtr; return temp; }

	//	ptrdiff_t operator-(const triangle_iterator& rawIterator){ return std::distance(rawIterator.getPtr(), this->getPtr()); }

	//	geometry::triangle<float, 3>& operator*(){ return *m_ptr; }
	//	const geometry::triangle<float, 3>& operator*()const{ return *m_ptr; }
	//	geometry::triangle<float, 3>* operator->(){ return m_ptr; }

	//	geometry::triangle<float, 3>* getPtr()const{ return m_ptr; }
	//	const geometry::triangle<float, 3>* getConstPtr()const{ return m_ptr; }

	//protected:

	//	int row;
	//	int col;
	//	int min_row;
	//	int max_row;
	//	int min_col; 
	//	int max_col;
	//	mesh3f_region* region;

	//	geometry::triangle<float, 3>* m_ptr;
	//};









}


#endif _GEOMETRY_MESH_H_
