#include "bspline_utils.h"


inline bspline::triangle3f bspline::mesh3f::get_triangle(int index)
{
	return bspline::triangle3f(points[elements[index][0]], points[elements[index][1]], points[elements[index][2]]);
}

VEC3F bspline::mesh3f::get_vertex(int element, int v){
	return points[elements[element][v]];
}


std::vector<bspline::triangle3f>* bspline::mesh3f::triangles()
{
	std::vector<bspline::triangle3f>* result = new std::vector<bspline::triangle3f>();
	for (int i = 0; i < elements.size(); i++)
	{
		result->push_back(get_triangle(i));
	}
	return result;
}

void bspline::mesh3f::get_mesh_range(VEC3F* min, VEC3F* max){
	*min = points[0];
	*max = points[0];
	for (int p = 1; p<points.size(); p++){
		min->x = fminf(min->x, points[p].x);
		min->y = fminf(min->y, points[p].y);
		min->z = fminf(min->z, points[p].z);
		max->x = fmaxf(max->x, points[p].x);
		max->y = fmaxf(max->y, points[p].y);
		max->z = fmaxf(max->z, points[p].z);
	}
}

inline void bspline::mesh3f::get_element_range(int i, VEC3F* min, VEC3F* max){
	*min = points[elements[i][0]];
	*max = points[elements[i][0]];
	for (int v = 1; v<3; v++){
		min->x = fminf(min->x, points[elements[i][v]].x);
		min->y = fminf(min->y, points[elements[i][v]].y);
		min->z = fminf(min->z, points[elements[i][v]].z);
		max->x = fmaxf(max->x, points[elements[i][v]].x);
		max->y = fmaxf(max->y, points[elements[i][v]].y);
		max->z = fmaxf(max->z, points[elements[i][v]].z);
	}
}

void bspline::mesh3f::build_nodes(int nx, int ny){

	nodes = lattice<std::vector<int>>(nx, ny);
	node_x = std::vector<float>(nx + 1);
	node_y = std::vector<float>(ny + 1);

	VEC3F min, max;
	get_mesh_range(&min, &max);

	for (int i = 0; i <= nx; i++){
		node_x[i] = (1 - (float)i / nx) * min.x + (float)i / nx * max.x;
	}
	for (int i = 0; i <= ny; i++){
		node_y[i] = (1 - (float)i / ny) * min.y + (float)i / ny * max.y;
	}

	for (int e = 0; e < elements.size(); e++){
		get_element_range(e, &min, &max);

		for (int i = 0; i < nx; i++){
			if (min.x <= node_x[i + 1] && max.x >= node_x[i])
			{
				for (int j = 0; j < ny; j++){
					if (min.y <= node_y[j + 1] && max.y >= node_y[j])
					{
						nodes.GetPoint(i, j).push_back(e);
					}
				}
			}
		}
	}
}

