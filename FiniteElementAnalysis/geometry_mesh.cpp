
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
	this->nodes_triangles = lattice<std::vector<int>>(nx, ny);
	this->nodes_edges = lattice<std::vector<int>>(nx, ny);
	this->nodes_points = lattice<std::vector<int>>(nx, ny);

	for (int e = 0; e < this->mesh->elements.size(); e++)
	{
		this->all_triangles.push_back(mesh_triangle<float>(mesh->get_vertex(e, 0), mesh->get_vertex(e, 1), mesh->get_vertex(e, 2)));
		this->all_edges.push_back(mesh_edge<float>(mesh->get_vertex(e, 0), mesh->get_vertex(e, 1)));
		this->all_edges.push_back(mesh_edge<float>(mesh->get_vertex(e, 1), mesh->get_vertex(e, 2)));
		this->all_edges.push_back(mesh_edge<float>(mesh->get_vertex(e, 2), mesh->get_vertex(e, 0)));
	}


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

	float t_min_x, t_min_y, t_max_x, t_max_y;


	//triangles
	for (int e = 0; e < this->all_triangles.size(); e++){
		mesh_triangle<float> t = all_triangles[e];

		t_min_x = (t.p1[0] < t.p2[0]) ? t.p1[0] : t.p2[0];
		if (t.p3[0] < t_min_x) t_min_x = t.p3[0];

		t_min_y = (t.p1[1] < t.p2[1]) ? t.p1[1] : t.p2[1];
		if (t.p3[1] < t_min_y) t_min_y = t.p3[1];

		t_max_x = (t.p1[0] < t.p2[0]) ? t.p2[0] : t.p1[0];
		if (t.p3[0] > t_min_x) t_min_x = t.p3[0];
		
		t_max_y = (t.p1[1] < t.p2[1]) ? t.p2[1] : t.p1[1];
		if (t.p3[1] > t_min_y) t_min_y = t.p3[1];

		for (int i = 0; i < nx; i++){
			if (t_min_x <= node_x[i + 1] && t_max_x >= node_x[i])
			{
				for (int j = 0; j < ny; j++){
					if (t_min_y <= node_y[j + 1] && t_max_y >= node_y[j])
					{
						nodes_triangles.GetPoint(i, j).push_back(e);
					}
				}
			}
		}
	}

	//edges
	for (int e = 0; e < this->all_edges.size(); e++){
		mesh_edge<float> t = all_edges[e];

		t_min_x = (t.p1[0] < t.p2[0]) ? t.p1[0] : t.p2[0];
		t_min_y = (t.p1[1] < t.p2[1]) ? t.p1[1] : t.p2[1];

		t_max_x = (t.p1[0] < t.p2[0]) ? t.p2[0] : t.p1[0];
		t_max_y = (t.p1[1] < t.p2[1]) ? t.p2[1] : t.p1[1];

		for (int i = 0; i < nx; i++){
			if (t_min_x <= node_x[i + 1] && t_max_x >= node_x[i])
			{
				for (int j = 0; j < ny; j++){
					if (t_min_y <= node_y[j + 1] && t_max_y >= node_y[j])
					{
						nodes_edges.GetPoint(i, j).push_back(e);
					}
				}
			}
		}
	}

	//points
	for (int e = 0; e < this->mesh->points.size(); e++){
		vector<float, 3> p = this->mesh->points[e];
		for (int i = 0; i < nx; i++){
			if (p[0] <= node_x[i + 1] && p[0] >= node_x[i])
			{
				for (int j = 0; j < ny; j++){
					if (p[1] <= node_y[j + 1] && p[1] >= node_y[j])
					{
						nodes_points.GetPoint(i, j).push_back(e);
					}
				}
			}
		}
	}
}



void geometry::mesh3f_region::set_region(float low_x, float high_x, float low_y, float high_y)
{

	int min_i = (int)(nx * (low_x - min_x) / (max_x - min_x));
	int min_j = (int)(ny * (low_y - min_y) / (max_y - min_y));
	int max_i = (int)(nx * (high_x - min_x) / (max_x - min_x));
	int max_j = (int)(ny * (high_y - min_y) / (max_y - min_y));

	if (min_i < 0) min_i = 0;
	if (min_j < 0) min_j = 0;
	if (max_i >= nodes_triangles.rows) max_i = nodes_triangles.rows - 1;
	if (max_j >= nodes_triangles.cols) max_j = nodes_triangles.cols - 1;

	if (min_i != min_x_index || min_j != min_y_index || max_i != max_x_index || max_j != max_y_index)
	{
		min_x_index = min_i;
		min_y_index = min_j;
		max_x_index = max_i;
		max_y_index = max_j;

		distinct_triangles.clear();
		distinct_edges.clear();
		distinct_points.clear();

		for (int i = min_x_index; i <= max_x_index; i++)
		{
			for (int j = min_y_index; j <= max_y_index; j++)
			{
				distinct_triangles.insert(distinct_triangles.end(), nodes_triangles.GetPoint(i, j).begin(), nodes_triangles.GetPoint(i, j).end());
				distinct_edges.insert(distinct_edges.end(), nodes_edges.GetPoint(i, j).begin(), nodes_edges.GetPoint(i, j).end());
				distinct_points.insert(distinct_points.end(), nodes_points.GetPoint(i, j).begin(), nodes_points.GetPoint(i, j).end());
			}
		}
		std::sort(distinct_triangles.begin(), distinct_triangles.end());
		auto last_triangle = std::unique(distinct_triangles.begin(), distinct_triangles.end());
		distinct_triangles.erase(last_triangle, distinct_triangles.end());

		std::sort(distinct_edges.begin(), distinct_edges.end());
		auto last_edge = std::unique(distinct_edges.begin(), distinct_edges.end());
		distinct_edges.erase(last_edge, distinct_edges.end());

		std::sort(distinct_points.begin(), distinct_points.end());
		auto last_point = std::unique(distinct_points.begin(), distinct_points.end());
		distinct_points.erase(last_point, distinct_points.end());
	}
}
