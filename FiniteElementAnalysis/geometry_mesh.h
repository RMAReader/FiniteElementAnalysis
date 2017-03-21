#ifndef _GEOMETRY_MESH_H_
#define _GEOMETRY_MESH_H_

#include <vector>
#include "geometry_vector.h"
#include "geometry_lattice.h"
#include "geometry_bspline_surface.h"


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

		mesh3f(geometry::bspline::surface<geometry::vector<float, 3 >,double, double>& surface, int nx, int ny);


		//inline triangle3f get_triangle(int index);
		geometry::vector<float, 3> get_vertex(int element, int v);
		//std::vector<triangle3f>* triangles();

		void get_mesh_range(geometry::vector<float, 3>& min, geometry::vector<float, 3>& max);
		inline void get_element_range(int i, geometry::vector<float, 3>& min, geometry::vector<float, 3>& max);
		void build_nodes(int nx, int ny);
		void build_grid_z(float, float);
	};


}


#endif _GEOMETRY_MESH_H_
