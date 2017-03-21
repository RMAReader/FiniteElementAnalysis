
#include "geometry.h"


//inline bspline::triangle3f bspline::mesh3f::get_triangle(int index)
//{
//	return bspline::triangle3f(points[elements[index][0]], points[elements[index][1]], points[elements[index][2]]);
//}

geoVEC3F geometry::mesh3f::get_vertex(int element, int v){
	return points[elements[element][v]];
}


//std::vector<bspline::triangle3f>* bspline::mesh3f::triangles()
//{
//	std::vector<bspline::triangle3f>* result = new std::vector<bspline::triangle3f>();
//	for (int i = 0; i < elements.size(); i++)
//	{
//		result->push_back(get_triangle(i));
//	}
//	return result;
//}

void geometry::mesh3f::get_mesh_range(geoVEC3F& min, geoVEC3F& max){
	min = points[0];
	max = points[0];
	for (int p = 1; p < points.size(); p++){
		for (int i = 0; i < 3; i++)
		{
			min[i] = fminf(min[i], points[p][i]);
			max[i] = fmaxf(max[i], points[p][i]);
		}
	}
}

inline void geometry::mesh3f::get_element_range(int e, geoVEC3F& min, geoVEC3F& max){
	min = points[elements[e][0]];
	max = points[elements[e][0]];
	for (int v = 1; v < 3; v++){
		for (int i = 0; i < 3; i++)
		{
			min[i] = fminf(min[i], points[elements[e][v]][i]);
			max[i] = fmaxf(max[i], points[elements[e][v]][i]);
		}
	}
}

void geometry::mesh3f::build_nodes(int nx, int ny){

	nodes = lattice<std::vector<int>>(nx, ny);
	node_x = std::vector<float>(nx + 1);
	node_y = std::vector<float>(ny + 1);

	geoVEC3F min, max;
	get_mesh_range(min, max);

	for (int i = 0; i <= nx; i++){
		node_x[i] = (1 - (float)i / nx) * min[0] + (float)i / nx * max[0];
	}
	for (int i = 0; i <= ny; i++){
		node_y[i] = (1 - (float)i / ny) * min[1] + (float)i / ny * max[1];
	}

	for (int e = 0; e < elements.size(); e++){
		get_element_range(e, min, max);

		for (int i = 0; i < nx; i++){
			if (min[0] <= node_x[i + 1] && max[0] >= node_x[i])
			{
				for (int j = 0; j < ny; j++){
					if (min[1] <= node_y[j + 1] && max[1] >= node_y[j])
					{
						nodes.GetPoint(i, j).push_back(e);
					}
				}
			}
		}
	}
}


geometry::mesh3f::mesh3f(geoSURFACE3F& surface, int nx, int ny)
{

	geoCURVE3F curve_y;

	float u, v;
	for (int i = 0; i < nx; i++){

		u = (1 - (float)i / (nx - 1)) * surface._knotx.minParam() + (float)i / (nx - 1) * surface._knotx.maxParam();
		curve_y = surface.curve_y(u);

		for (int j = 0; j < ny; j++){
			v = (1 - (float)j / (ny - 1)) * surface._knoty.minParam() + (float)j / (ny - 1) * surface._knoty.maxParam();

			points.push_back(curve_y.evaluate(v));

		}
	}

	int p1, p2, p3, p4, count;
	elements.reserve(2 * nx * ny);
	for (int i = ny; i < nx * ny; i += ny){
		for (int j = 1; j < ny; j++){
			p1 = i + j;
			p2 = p1 - 1;
			p3 = p1 - ny;
			p4 = p3 - 1;
			elements.push_back({ p1, p2, p3 });
			elements.push_back({ p4, p3, p2 });
		}
	}
}