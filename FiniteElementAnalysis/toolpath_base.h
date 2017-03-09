

#ifndef _TOOLPATH_H_
#define _TOOLPATH_H_



#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <time.h>
#include <iomanip>

#include "bspline_utils.h"
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
	std::vector<VEC3F> points;

	virtual void calculate();
	void save_gcode();

	void bore_hole(circle2f& c, float safe_z, float max_z, float min_z, float step_z);
	
	void cut_curve_add_point(std::vector<VEC3F>* path, VEC3F new_p, VEC2F& tangent, float resolution);
	void cut_polyline_const_z(std::vector<VEC3F>* path, CURVE2F& c, float z, double minParam, double maxParam);
	void cut_curve_const_z(std::vector<VEC3F>* path, float tool_diameter, CURVE2F& c, float z, double minParam, double maxParam, float resolution);
	void cut_curve_const_z(std::vector<VEC3F>* path, float tool_diameter, bspline::offsetCurve2f& c, float z, double minParam, double maxParam, float resolution);
	void offset_curve(std::vector<VEC3F>* path, float tool_diameter, bool deep_corners);
	void join_curves_external(std::vector<VEC3F>& path1, std::vector<VEC3F>& path2);
	void close_curve_external(std::vector<VEC3F>& path);

	void trim_corner_block(CURVE2F* c_bout, CURVE2F* bout, CURVE2F* block, float block_top_z, float block_bottom_z, float step, float resolution);
	void trim_end_block(CURVE2F* rib, CURVE2F* centre_line, CURVE2F* block, float block_top_z, float block_bottom_z, float step, float resolution);

	void rough_surface_scanning_stl(std::vector<VEC3F>* path, float tool_diameter, float step_x, float step_y, float minimum_z, bspline::mesh3f*);
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


#endif _TOOLPATH_H_