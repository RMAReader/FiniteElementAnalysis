
#include "toolpath_base.h"

void toolpath_TrimBlocksCentreBout::calculate()
{
	using namespace bspline;

	float safe_z = parameters["safe_z"];
	float block_top_z = parameters["block_top_z"];
	float block_bottom_z = parameters["block_bottom_z"];
	float step_z = parameters["step_z"];
	float resolution = 0.1;

	std::unordered_map<std::string, CURVE2F*>& c = violin->ribs.curves;
	points.push_back(VEC3F(0, 0, safe_z));

	trim_corner_block(c["rib_internal_centre_bout_trbl"], c["rib_internal_lower_bout"], c["block_trbl_lower_bout"], block_top_z, block_bottom_z, step_z, resolution);
	trim_corner_block(c["rib_internal_centre_bout_trbl"], c["rib_internal_upper_bout"], c["block_trbl_upper_bout"], block_top_z, block_bottom_z, step_z, resolution);
	trim_corner_block(c["rib_internal_centre_bout_bass"], c["rib_internal_upper_bout"], c["block_bass_upper_bout"], block_top_z, block_bottom_z, step_z, resolution);
	trim_corner_block(c["rib_internal_centre_bout_bass"], c["rib_internal_lower_bout"], c["block_bass_lower_bout"], block_top_z, block_bottom_z, step_z, resolution);

}

