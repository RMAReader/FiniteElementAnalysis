#include "toolpath_base.h"

void toolpath_BackRough::calculate()
{
	using namespace bspline;

	float safe_z = parameters["safe_z"];
	float minimum_z = parameters["minimum_z"];
	float step_x = parameters["step_x"];
	float step_y = parameters["step_y"];
	float tool_diameter = parameters["tool_diameter"];
	int nx = (int)parameters["mesh_nx"];
	int ny = (int)parameters["mesh_ny"];

	SURFACE3F* s = violin->back->surfaces["exterior"];
	mesh3f* mesh = bspline::build_mesh(s, nx, ny);

	points.push_back(VEC3F(0, 0, safe_z));

	rough_surface_scanning_stl(&points, tool_diameter, step_x, step_y, minimum_z, mesh);
}