#include "toolpath_base.h"

void toolpath_RibMould::calculate()
{
	//zero z is set to top of rib mould base

	float safe_z = parameters["safe_z"];
	float mould_top_z = parameters["mould_top_z"];
	float mould_bottom_z = parameters["mould_bottom_z"];
	float step_z = parameters["step_z"];


	//bore location holes
	for (int i = 0; i < violin->ribs.rib_mould_locator_holes.size(); i++)
	{
		geoVEC2F c;
		c[0] = violin->ribs.rib_mould_locator_holes[i].centre.x;
		c[1] = violin->ribs.rib_mould_locator_holes[i].centre.y;
		float d = violin->ribs.rib_mould_locator_holes[i].diameter;

		bore_hole(c, d, safe_z, mould_top_z, mould_bottom_z, step_z);
	}

	float resolution = 1.0;


	//cut out central clamping hole
	std::vector<geoVEC3F> clamp_hole_layer;
	std::unordered_map<std::string, geoCURVE2F>& c = violin->ribs.curves;

	cut_curve_const_z(&clamp_hole_layer, tool->diameter, c["clamp_hole"], 0, c["clamp_hole"].minParam(), c["clamp_hole"].maxParam(), resolution);

	points.push_back(geoVEC3F(std::array < float, 3 > {{clamp_hole_layer[0][0], clamp_hole_layer[0][1], safe_z}}));
	float z = mould_top_z;
	while (z > mould_bottom_z)
	{
		z -= step_z;
		if (z < mould_bottom_z){ z = mould_bottom_z; }

		for (int i = 0; i < clamp_hole_layer.size(); i++){
			points.push_back(geoVEC3F(std::array < float, 3 > {{clamp_hole_layer[i][0], clamp_hole_layer[i][1], z}}));
		}
	}
	int n = clamp_hole_layer.size() - 1;
	points.push_back(geoVEC3F(std::array < float, 3 > {{clamp_hole_layer[n][0], clamp_hole_layer[n][1], safe_z}}));



	//cut outline - clockwise from lower bout treble corner
	float error = 1E-4;

	//lower bout
	bool intersect1, intersect2, intersect3, intersect4, intersect5, intersect6;
	bool intersect7, intersect8, intersect9, intersect10, intersect11, intersect12;
	geoVEC2D param1, param2, param3, param4;

	intersect1 = IntersectParam(c["rib_internal_lower_bout"], c["block_trbl_lower_bout"], error, param1, false);
	intersect2 = IntersectParam(c["rib_internal_lower_bout"], c["centre_line"], error, param2, false);
	intersect4 = IntersectParam(c["rib_internal_lower_bout"], c["block_bass_lower_bout"], error, param4, false);
	if (intersect1 && intersect2 && intersect4){
		geoCURVE2F rib_lower_trbl = c["rib_internal_lower_bout"]; rib_lower_trbl.trim_curve(param1[0], param2[0]);
		geoCURVE2F rib_lower_bass = c["rib_internal_lower_bout"]; rib_lower_bass.trim_curve(param2[0], param4[0]);

		intersect2 = IntersectParam(rib_lower_trbl, c["block_bottom"], error, param2, false);
		intersect3 = IntersectParam(rib_lower_bass, c["block_bottom"], error, param3, false);
	}
	//centre bout bass
	geoVEC2D param5, param6;
	intersect5 = IntersectParam(c["rib_internal_centre_bout_bass"], c["block_bass_lower_bout"], error, param5, false);
	intersect6 = IntersectParam(c["rib_internal_centre_bout_bass"], c["block_bass_upper_bout"], error, param6, false);

	//upper bout
	geoVEC2D param7, param8, param9, param10;

	intersect7 = IntersectParam(c["rib_internal_upper_bout"], c["block_bass_upper_bout"], error, param7, false);
	intersect8 = IntersectParam(c["rib_internal_upper_bout"], c["centre_line"], error, param8, false);
	intersect10 = IntersectParam(c["rib_internal_upper_bout"], c["block_trbl_upper_bout"], error, param10, false);
	if (intersect7 && intersect8 && intersect10){
		geoCURVE2F rib_upper_bass = c["rib_internal_upper_bout"]; rib_upper_bass.trim_curve(param7[0], param8[0]);
		geoCURVE2F rib_upper_trbl = c["rib_internal_upper_bout"]; rib_upper_trbl.trim_curve(param8[0], param10[0]);

		intersect8 = IntersectParam(rib_upper_bass, c["block_neck"], error, param8, false);
		intersect9 = IntersectParam(rib_upper_trbl, c["block_neck"], error, param9, false);
	}

	//centre bout treble
	geoVEC2D param11, param12;
	intersect11 = IntersectParam(c["rib_internal_centre_bout_trbl"], c["block_trbl_upper_bout"], error, param11, false);
	intersect12 = IntersectParam(c["rib_internal_centre_bout_trbl"], c["block_trbl_lower_bout"], error, param12, false);

	if (intersect1 && intersect2 && intersect3 && intersect4 && intersect5 && intersect6 &&
		intersect7 && intersect8 && intersect9 && intersect10 && intersect11 && intersect12)
	{
		std::vector<geoVEC3F> layer1, layer2, layer3, layer4, layer5, layer6, layer7, layer8, layer9, layer10, layer11, layer12;

		cut_curve_const_z(&layer1, 0, c["rib_internal_lower_bout"], 0, param1[0], param2[0], resolution);
		cut_polyline_const_z(&layer2, c["block_bottom"], 0, param2[1], param3[1]);
		cut_curve_const_z(&layer3, 0, c["rib_internal_lower_bout"], 0, param3[0], param4[0], resolution);

		cut_polyline_const_z(&layer4, c["block_bass_lower_bout"], 0, param4[1], param5[1]);
		cut_curve_const_z(&layer5, 0, c["rib_internal_centre_bout_bass"], 0, param5[0], param6[0], resolution);
		cut_polyline_const_z(&layer6, c["block_bass_upper_bout"], 0, param6[1], param7[1]);

		cut_curve_const_z(&layer7, 0, c["rib_internal_upper_bout"], 0, param7[0], param8[0], resolution);
		cut_polyline_const_z(&layer8, c["block_neck"], 0, param8[1], param9[1]);
		cut_curve_const_z(&layer9, 0, c["rib_internal_upper_bout"], 0, param9[0], param10[0], resolution);

		cut_polyline_const_z(&layer10, c["block_trbl_upper_bout"], 0, param10[1], param11[1]);
		cut_curve_const_z(&layer11, 0, c["rib_internal_centre_bout_trbl"], 0, param11[0], param12[0], resolution);
		cut_polyline_const_z(&layer12, c["block_trbl_lower_bout"], 0, param12[1], param1[1]);

		//rib outlines
		bool deep_corners = false;
		offset_curve(&layer1, tool->diameter, deep_corners);
		offset_curve(&layer3, tool->diameter, deep_corners);
		offset_curve(&layer5, tool->diameter, deep_corners);
		offset_curve(&layer7, tool->diameter, deep_corners);
		offset_curve(&layer9, tool->diameter, deep_corners);
		offset_curve(&layer11, tool->diameter, deep_corners);

		//block pockets
		deep_corners = true;
		offset_curve(&layer2, tool->diameter, deep_corners);
		offset_curve(&layer4, tool->diameter, deep_corners);
		offset_curve(&layer6, tool->diameter, deep_corners);
		offset_curve(&layer8, tool->diameter, deep_corners);
		offset_curve(&layer10, tool->diameter, deep_corners);
		offset_curve(&layer12, tool->diameter, deep_corners);

		std::vector<geoVEC3F> layer;
		join_curves_external(layer, layer1);
		join_curves_external(layer, layer2);
		join_curves_external(layer, layer3);
		join_curves_external(layer, layer4);
		join_curves_external(layer, layer5);
		join_curves_external(layer, layer6);
		join_curves_external(layer, layer7);
		join_curves_external(layer, layer8);
		join_curves_external(layer, layer9);
		join_curves_external(layer, layer10);
		join_curves_external(layer, layer11);
		join_curves_external(layer, layer12);
		close_curve_external(layer);

		points.push_back(geoVEC3F(std::array < float, 3 > {{layer.front()[0], layer.front()[1], safe_z}}));
		z = mould_top_z;
		while (z > mould_bottom_z)
		{
			z -= step_z;
			if (z < mould_bottom_z){ z = mould_bottom_z; }

			for (int i = 0; i < layer.size(); i++){
				points.push_back(geoVEC3F(std::array < float, 3 > {{layer[i][0], layer[i][1], z}}));
			}
		}
		points.push_back(geoVEC3F(std::array < float, 3 > {{layer.back()[0], layer.back()[1], safe_z}}));
	}
}
