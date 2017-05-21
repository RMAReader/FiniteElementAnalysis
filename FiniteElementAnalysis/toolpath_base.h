

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



/*
interface which provides function to evaluate z coordinate of toolpath, given x and y coordinates
*/
class toolpath_height_base
{
public:
	virtual float get_height(float x, float y);
	virtual float default_height();
};

class toolpath_mesh_height : public toolpath_height_base
{
private:
	geometry::mesh3f_region region;
	float r;
	float default_z;
public:

	toolpath_mesh_height(geometry::mesh3f* mesh, int nx, int ny, float r, float default_z);
	float get_height(float x, float y);
	float default_height();
};


template <class T>
struct scan_line
{
	T y;
	std::vector<T> x;
	scan_line();
	scan_line(T y){ this->y = y; }
	scan_line(T y, std::vector<T>& x){ this->y = y; this->x = x; }
};



template <class T, int N>
class range
{
private:
	std::vector<geometry::vector<T, N>> _range;
	geometry::vector<T, N> min;
	geometry::vector<T, N> max;
	T step;
	bool include_back;

	void build_range()
	{
		if (include_back)
		{
			int n = (int)((max - min).L2norm()) / step + 1;
			_range.resize(n + 1);
			for (int i = 0; i <= n; i++)
			{
				_range[i] = (1 - (T)i / n) * min + (T)i / n * max;
			}
		}
		else{
			int n = (int)((max - min).L2norm()) / step + 1;
			_range.resize(n);
			for (int i = 0; i < n; i++)
			{
				_range[i] = (1 - (T)i / n) * min + (T)i / n * max;
			}
		}
	}

public:

	range(geometry::vector<T, N> min, geometry::vector<T, N> max, float step, bool include_back)
	{
		this->min = min;
		this->max = max;
		this->step = step;
		this->include_back = include_back;
		build_range();
	}
	
	//std::vector<geometry::vector<T, N>>* get_range(){return &_range;}
	typedef typename std::vector<geometry::vector<T, N>>::iterator range_iterator;

	range_iterator begin(){ return _range.begin(); }
	range_iterator end(){ return _range.end(); }
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

	static void build_scanlines(std::vector<scan_line<float>>&, std::vector<geometry::line<float, 2>>& border, float step_y);
	static void scanline_path_2D(std::vector<std::vector<geoVEC2F>>* path, std::vector<scan_line<float>> scan_lines, float tool_diameter);

	static void scanning_path_3D(std::vector<geoVEC3F>* path, std::vector<std::vector<geoVEC2F>>& path_2D, float step, toolpath_height_base* h);

	static std::vector<float> toolpath_base::get_range(float min, float max, int n);
	static std::vector<float> toolpath_base::get_range(float min, float max, float step);

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