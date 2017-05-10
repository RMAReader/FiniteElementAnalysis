
#include "geometry.h"


geometry::mesh3f::mesh3f()
{}

geoVEC3F geometry::mesh3f::get_vertex(int element, int v){
	return points[elements[element][v]];
}



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



void geometry::mesh3f::add_triangle(triangle<float, 3>& t)
{
	int i = points.size();
	points.push_back(t.p1);
	points.push_back(t.p2);
	points.push_back(t.p3);
	elements.push_back({i, i+1, i+2});
}



geometry::mesh3f_region::mesh3f_region()
{
}

geometry::mesh3f_region::mesh3f_region(mesh3f* mesh, int nx, int ny)
{
	this->build_nodes(mesh, nx, ny);
}

void geometry::mesh3f_region::build_nodes(mesh3f* mesh, int nx, int ny)
{
	this->mesh = mesh;
	this->nodes = lattice<std::vector<int>>(nx, ny);

	geoVEC3F min, max;
	this->mesh->get_mesh_range(min, max);
	this->max_x = max[0];
	this->max_y = max[1];
	this->min_x = min[0];
	this->min_y = min[1];
	this->nx = nx;
	this->ny = ny;

	std::vector<float> node_x(nx + 1);
	std::vector<float> node_y(ny + 1);

	for (int i = 0; i <= nx; i++){
		node_x[i] = (1 - (float)i / nx) * min[0] + (float)i / nx * max[0];
	}
	for (int i = 0; i <= ny; i++){
		node_y[i] = (1 - (float)i / ny) * min[1] + (float)i / ny * max[1];
	}

	for (int e = 0; e < this->mesh->elements.size(); e++){
		this->mesh->get_element_range(e, min, max);

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



void geometry::mesh3f_region::set_region(float low_x, float high_x, float low_y, float high_y)
{
	this->points.clear();
	this->edges.clear();
	this->triangles.clear();
	
	int min_x_index = (int)(nx * (low_x - min_x) / (max_x - min_x));
	int min_y_index = (int)(ny * (low_y - min_y) / (max_y - min_y));
	int max_x_index = (int)(nx * (high_x - min_x) / (max_x - min_x));
	int max_y_index = (int)(ny * (high_y - min_y) / (max_y - min_y));

	if (min_x_index < 0) min_x_index = 0;
	if (min_y_index < 0) min_y_index = 0;
	if (max_x_index >= nodes.rows) max_x_index = nodes.rows-1;
	if (max_y_index >= nodes.cols) max_y_index = nodes.cols-1;

	for (int i = min_x_index; i <= max_x_index; i++)
	{
		for (int j = min_y_index; j <= max_y_index; j++)
		{
			for each(auto e in nodes.GetPoint(i, j))
			{
				this->triangles.push_back(geometry::triangle<float, 3>(mesh->get_vertex(e, 0), mesh->get_vertex(e, 1), mesh->get_vertex(e, 2)));
				this->edges.push_back(geometry::line<float, 3>(mesh->get_vertex(e, 0), mesh->get_vertex(e, 1)));
				this->edges.push_back(geometry::line<float, 3>(mesh->get_vertex(e, 1), mesh->get_vertex(e, 2)));
				this->edges.push_back(geometry::line<float, 3>(mesh->get_vertex(e, 2), mesh->get_vertex(e, 0)));
				this->points.push_back(geometry::vector<float, 3>(mesh->get_vertex(e, 0)));
				this->points.push_back(geometry::vector<float, 3>(mesh->get_vertex(e, 1)));
				this->points.push_back(geometry::vector<float, 3>(mesh->get_vertex(e, 2)));
			}
		}
	}
}
