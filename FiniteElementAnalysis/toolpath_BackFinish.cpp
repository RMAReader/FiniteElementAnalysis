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
	geoCURVE2F b1, b2, b3, b4;
	geometry::bspline::Project(b1, s.curve_x(s._knoty.minParam()));
	geometry::bspline::Project(b2, s.curve_y(s._knotx.minParam()));
	geometry::bspline::Project(b3, s.curve_x(s._knoty.maxParam()));
	geometry::bspline::Project(b4, s.curve_y(s._knotx.maxParam()));

	auto c1 = geometry::bspline::ToPolyLine(b1, step_y);
	auto c2 = geometry::bspline::ToPolyLine(b2, step_y);
	auto c3 = geometry::bspline::ToPolyLine(b3, step_y);
	auto c4 = geometry::bspline::ToPolyLine(b4, step_y);

	std::vector < geometry::line<float, 2>> border;
	border.insert(border.end(), c1.begin(), c1.end());
	border.insert(border.end(), c2.begin(), c2.end());
	border.insert(border.end(), c3.begin(), c3.end());
	border.insert(border.end(), c4.begin(), c4.end());



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