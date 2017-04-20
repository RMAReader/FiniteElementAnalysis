

#ifndef _TOOLPATH_H_
#define _TOOLPATH_H_



#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <time.h>
#include <iomanip>

//#include "bspline_utils.h"
#include "geometry.h"
#include "violin_model.h"


class cnc_tool
{
public:
	std::string name;
	std::string type;
	float diameter;
};


class toolpath_base
{
public:
	toolpath_base();
	~toolpath_base();

	violin_model* violin;
	cnc_tool* tool;
	std::string path_name;
	std::string gcode_filepath;
	std::unordered_map<std::string, float> parameters;
	std::vector<geoVEC3F> points;

	virtual void calculate();
	void save_gcode();

	void bore_hole(geoVEC2F& c, float d, float safe_z, float max_z, float min_z, float step_z);
	
	void cut_curve_add_point(std::vector<geoVEC3F>* path, geoVEC3F new_p, geoVEC2F& tangent, float resolution);
	void cut_polyline_const_z(std::vector<geoVEC3F>* path, geoCURVE2F& c, float z, double minParam, double maxParam);
	//void cut_curve_const_z(std::vector<geoVEC3F>* path, float tool_diameter, geoCURVE2F& c, float z, double minParam, double maxParam, float resolution);
	void cut_curve_const_z(std::vector<geoVEC3F>* path, float tool_diameter, geoCURVE2F& c, float z, double minParam, double maxParam, float resolution);
	void offset_curve(std::vector<geoVEC3F>* path, float tool_diameter, bool deep_corners);
	void join_curves_external(std::vector<geoVEC3F>& path1, std::vector<geoVEC3F>& path2);
	void close_curve_external(std::vector<geoVEC3F>& path);

	void trim_corner_block(geoCURVE2F& c_bout, geoCURVE2F& bout, geoCURVE2F& block, float block_top_z, float block_bottom_z, float step_z, float resolution);
	void trim_end_block(geoCURVE2F& rib, geoCURVE2F& centre_line, geoCURVE2F& block, float block_top_z, float block_bottom_z, float step_z, float resolution);

	void finish_surface_point_cloud(std::vector<geoVEC3F>* path, float tool_diameter, float step_x, float step_y, float minimum_z, geometry::mesh3f& mesh);
	void finish_surface_scanning_stl(std::vector<geoVEC3F>* path, float tool_diameter, float step_x, float step_y, float minimum_z, geometry::mesh3f&);
	void rough_surface_grid(std::vector<geoVEC3F>* path, float tool_diameter, float step_x, float step_y, float minimum_z, float safe_z, float z_margin, geometry::mesh3f&);
};

class toolpath_RibMouldBaseJig : public toolpath_base
{
public:
	void calculate();
};

class toolpath_RibMould : public toolpath_base
{
public:
	void calculate();
};

class toolpath_TrimBlocksCentreBout : public toolpath_base
{
public:
	void calculate();
};


class toolpath_TrimBlocksEndBouts : public toolpath_base
{
public:
	void calculate();
};



class toolpath_BackRough : public toolpath_base
{
public:
	void calculate();
};

class toolpath_BackFinish : public toolpath_base
{
public:
	void calculate();
};


#endif _TOOLPATH_H_