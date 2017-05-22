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


	class mesh3f_region
	{
	private:
		mesh3f* mesh;
		
		float min_x;
		float min_y;
		float max_x;
		float max_y;
		int nx;
		int ny;

		std::vector<geometry::triangle<float, 3>> all_triangles;
		std::vector<geometry::line<float, 3>> all_edges;
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

		inline geometry::triangle<float, 3> get_triangle(int i){ return all_triangles[i]; }
		inline geometry::line<float, 3> get_edge(int i){ return all_edges[i]; }
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








	template <class T>
	struct mesh_edge
	{
		vector<T, 3> t;
		vector<T, 3> p1;
		T u1,u3;
		T a;
	public:
		mesh_edge(vector<T, 3>& v1, vector<T, 3>& v2)
		{
			p1 = v1;
			t = v2 - v1;

			u1 = (t[0] * t[0] + t[1] * t[1]) / t[2];

			a = t[0] * t[0] + t[1] * t[1] + u1 * u1;
		}

		bool min_height_sphere(T x, T y, T r, T& h)
		{
			vector<T, 2> p;  
			p[0] = p1[0] - x; 
			p[1] = p1[1] - y;

			T u2 = (p[0] * t[0] + p[1] * t[1]) / t[2];

			T b = 2 * (p[0] * t[0] + p[1] * t[1] + u1 * u2);
			T c = p[0] * p[0] + p[1] * p[1] + u2 * u2 - r * r;
			
			T d = b * b - 4 * a * c;
			T alpha;
			T h1 = 0;
			T h2 = 0;
			bool contact = false;
			if (d >= 0)
			{
				d = sqrtf(d);

				alpha = (-b + d) / (2 * a);
				if (0 <= alpha && alpha <= 1)
				{
					h1 = (alpha * u1 + u2) + p1[2] + alpha * t[2];
					contact = true;
				}
				alpha = (-b - d) / (2 * a);
				if (0 <= alpha && alpha <= 1)
				{
					h1 = (alpha * u1 + u2) + p1[2] + alpha * t[2];
					contact = true;
				}

			}
			if (contact){
				h = fmaxf(h1, h2);
				return true;
			}
			else{
				return false;
			}

		}

	};
}


#endif _GEOMETRY_MESH_H_
