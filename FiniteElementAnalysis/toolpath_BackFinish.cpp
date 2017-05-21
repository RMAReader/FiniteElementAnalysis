#include "toolpath_base.h"
#include <chrono>

void toolpath_BackFinish::calculate()
{

	if (violin == nullptr) { return; }

	if (violin->back.surfaces.find("exterior") == violin->back.surfaces.end()) return;

	float safe_z = parameters["safe_z"];
	float minimum_z = parameters["minimum_z"];
	float margin_z = parameters["margin_z"];
	float step_x = parameters["step_x"];
	float step_y = parameters["step_y"];
	int nx = (int)parameters["mesh_nx"];
	int ny = (int)parameters["mesh_ny"];

	geoSURFACE3F s = violin->back.surfaces["exterior"];
	geometry::mesh3f mesh(s, nx, ny);

	points.push_back(geoVEC3F(std::array < float, 3 > {{0, 0, safe_z}}));

	auto t1 = std::chrono::high_resolution_clock::now();

	//finish_surface_scanning_stl(&points, tool->diameter, step_x, step_y, minimum_z, mesh);
	

	//get border
	geoVEC3F min, max;
	mesh.get_mesh_range(min, max);
	geoVEC2F p1(std::array < float, 2 > {{min[0], min[1]}});
	geoVEC2F p2(std::array < float, 2 > {{min[0], max[1]}});
	geoVEC2F p3(std::array < float, 2 > {{max[0], max[1]}});
	geoVEC2F p4(std::array < float, 2 > {{max[0], min[1]}});
	std::vector < geometry::line<float, 2>> border
	{
		geometry::line<float, 2 >(p1, p2),
		geometry::line<float, 2 >(p2, p3),
		geometry::line<float, 2 >(p3, p4),
		geometry::line<float, 2 >(p4, p1),
	};

	//build scan lines
	std::vector<scan_line<float>> scan_lines;
	build_scanlines(scan_lines, border, step_y);

	//build 2d path
	std::vector<std::vector<geoVEC2F>> path_2D;
	scanline_path_2D(&path_2D, scan_lines, tool->diameter);

	//build 3d path
	toolpath_mesh_height h(&mesh, 70, 50, tool->diameter * 0.5, safe_z);
	scanning_path_3D(&points, path_2D, step_y, &h);

	auto t2 = std::chrono::high_resolution_clock::now();
	std::cout << "finish_surface_scanning_stl() took "
		<< std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count()
		<< " seconds\n";
}