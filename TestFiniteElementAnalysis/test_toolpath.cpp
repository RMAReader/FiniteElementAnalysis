#include "test_base.h"
#include "geometry.h"
#include "toolpath_base.h"

void toolpath_scanning_path_3D::run()
{
	name = "toolpath_scanning_path_3D";
	float tolerance = 0.000001;

	geometry::mesh3f mesh;
	geoVEC3F p1(std::array < float, 3 > {{0, 0, 0}});
	geoVEC3F p2(std::array < float, 3 > {{1, 0, 0}});
	geoVEC3F p3(std::array < float, 3 > {{0, 2, 0}});
	mesh.add_triangle(geometry::triangle<float, 3>(p1, p2, p3));
	
	int nx = 1;
	int ny = 1;
	float radius = 0.2;
	float default_z = 30;

	toolpath_mesh_height h(&mesh, nx, ny, radius, default_z);

	std::vector<std::vector<geoVEC2F>> path_2D;
	std::vector<geoVEC2F> path_section;
	path_section.push_back(geoVEC2F(std::array < float, 2 > {{-0.1, 0.1}}));
	path_section.push_back(geoVEC2F(std::array < float, 2 > {{0.0, 0.1}}));
	path_section.push_back(geoVEC2F(std::array < float, 2 > {{0.1, 0.1}}));
	path_section.push_back(geoVEC2F(std::array < float, 2 > {{0.2, 0.1}}));
	path_section.push_back(geoVEC2F(std::array < float, 2 > {{0.8, 0.1}}));
	path_section.push_back(geoVEC2F(std::array < float, 2 > {{0.9, 0.1}}));
	path_section.push_back(geoVEC2F(std::array < float, 2 > {{1.0, 0.1}}));
	path_section.push_back(geoVEC2F(std::array < float, 2 > {{1.1, 0.1}}));

	path_2D.push_back(path_section);

	std::vector<geoVEC3F> path_3D;
	toolpath_base::scanning_path_3D(&path_3D, path_2D, &h);
}