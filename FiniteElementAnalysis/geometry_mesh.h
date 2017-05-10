#ifndef _GEOMETRY_MESH_H_
#define _GEOMETRY_MESH_H_

#include <vector>
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
		
		lattice<std::vector<int>> nodes;
		float min_x;
		float min_y;
		float max_x;
		float max_y;
		int nx;
		int ny;

		std::vector<geometry::triangle<float, 3>> triangles;
		std::vector<geometry::line<float, 3>> edges;
		std::vector<geometry::vector<float, 3>> points;
	public:

		mesh3f_region();
		mesh3f_region(mesh3f* mesh, int nx, int ny);

		void build_nodes(mesh3f* mesh, int nx, int ny);
		void set_region(float min_x, float max_x, float min_y, float max_y);

		typedef std::vector<geometry::triangle<float, 3>>::iterator triangle_iterator;
		triangle_iterator begin_triangle() { return triangles.begin(); }
		triangle_iterator end_triangle()   { return triangles.end(); }

		typedef std::vector<geometry::line<float, 3>>::iterator edge_iterator;
		edge_iterator begin_edge() { return edges.begin(); }
		edge_iterator end_edge()   { return edges.end(); }

		typedef std::vector<geometry::vector<float, 3>>::iterator point_iterator;
		point_iterator begin_point() { return points.begin(); }
		point_iterator end_point()   { return points.end(); }

	};




	template <class T>
	struct mesh_edge
	{
		vector<T, 3> t;
		vector<T, 3> p1;
		T u1,u3;
		T a;
	public:
		void initialise(vector<T, 3>& v1, vector<T, 3>& v2)
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
