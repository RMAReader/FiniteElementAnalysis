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

	finish_surface_scanning_stl(&points, tool->diameter, step_x, step_y, minimum_z, mesh);
	
	auto t2 = std::chrono::high_resolution_clock::now();
	std::cout << "finish_surface_scanning_stl() took "
		<< std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count()
		<< " seconds\n";
}