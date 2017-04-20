
#include "toolpath_base.h"


toolpath_base::toolpath_base()
{
}


toolpath_base::~toolpath_base()
{
}

void toolpath_base::calculate()
{
}

void toolpath_base::save_gcode()
{
	using namespace std;
	fstream out;
	out.open(gcode_filepath, ios::out);
	if (out.is_open()){

		time_t rawtime;
		struct tm * timeinfo;
		char buffer[80];
		time(&rawtime);
		timeinfo = localtime(&rawtime);
		strftime(buffer, 80, "%c", timeinfo);

		out << "(" << buffer << ")" << endl;
		out << "(" << path_name << ")" << endl;
		out << "(" << tool->name << ")" << endl;

		out << "G90"; //Select absolute distance mode 
		out << "G49"; //Cancel tool length offset 
		out << "S" << fixed << setprecision(0) << parameters["spindle_speed"] << "M03" << endl;

		out << "G01";
		out << "Z" << fixed << setprecision(3) << parameters["safe_z"];
		out << "F" << fixed << setprecision(0) << parameters["z_feedrate"];
		out << endl;

		if (points.size() > 0)
		{

			out << "G01";
			out << "X" << fixed << setprecision(3) << points[0][0];
			out << "Y" << fixed << setprecision(3) << points[0][1];
			out << "F" << fixed << setprecision(0) << parameters["max_feedrate"];
			out << endl;

			for (int i = 0; i < points.size(); i++)
			{
				out << "G01";
				out << "X" << fixed << setprecision(3) << points[i][0];
				out << "Y" << fixed << setprecision(3) << points[i][1];
				out << "Z" << fixed << setprecision(3) << points[i][2];
				out << "F" << fixed << setprecision(0) << parameters["xy_feedrate"];
				out << endl;
			}
		}
		out << "M05"; //spindle off
		out << "M30"; //end program
		out << endl;

		out.close();
	}

}


void toolpath_base::bore_hole(geoVEC2F& c, float d, float safe_z, float max_z, float min_z, float step_z)
{
	if (tool == nullptr){ exit; }
	if (tool->diameter < 0){ exit; }


	float radius = (d - tool->diameter) / 2;
	if (radius > 0)
	{

		int n = 100;
		int j = 0;
		double PI = 3.141592653589793238463;
		float z = max_z;
		float theta = 0;

		//move to safe z and centre of circle
		points.push_back(geoVEC3F(std::array < float, 3 > {{c[0], c[1], safe_z}}));
		points.push_back(geoVEC3F(std::array < float, 3 > {{c[0], c[1], z}}));

		//spiral down to required depth
		while (z >= min_z)
		{
			theta = PI * 2 * j / n;
			z -= step_z / n;
			if (z >= min_z)
			{
				geoVEC3F p;
				p[0] = c[0] + radius * sinf(theta);
				p[1] = c[1] + radius * cosf(theta);
				p[2] = z;
				points.push_back(p);
				j++;
			}
		}

		//do one full circle at required depth
		for (int k = j; k <= j + n; k++)
		{
			theta = PI * 2 * k / n;
			geoVEC3F p;
			p[0] = c[0] + radius * sinf(theta);
			p[1] = c[1] + radius * cosf(theta);
			p[2] = min_z;
			points.push_back(p);
		}

		//move to centre of circle and safe z
		points.push_back(geoVEC3F(std::array < float, 3 > {{c[0], c[1], min_z}}));
		points.push_back(geoVEC3F(std::array < float, 3 > {{c[0], c[1], safe_z}}));
	}
}



void toolpath_base::cut_polyline_const_z(std::vector<geoVEC3F>* path, geoCURVE2F& c, float z, double minParam, double maxParam)
{
	double min_t, max_t;
	geoVEC2F new_p;

	if (minParam < maxParam){
		min_t = (c.minParam() < minParam) ? minParam : c.minParam();
		max_t = (c.maxParam() > maxParam) ? maxParam : c.maxParam();
		new_p = c.evaluate(min_t);
		path->push_back(geoVEC3F(std::array < float, 3 > {{new_p[0], new_p[1], z}}));
		for (int i = 0; i < c._knot.size(); i++){
			if (min_t < c._knot[i] && c._knot[i] < max_t){
				new_p = c.evaluate(c._knot[i]);
				path->push_back(geoVEC3F(std::array < float, 3 > {{new_p[0], new_p[1], z}}));
			}
		}
		new_p = c.evaluate(max_t);
		path->push_back(geoVEC3F(std::array < float, 3 > {{new_p[0], new_p[1], z}}));
	}
	else{
		max_t = (c.maxParam() > minParam) ? minParam : c.maxParam();
		min_t = (c.minParam() < maxParam) ? maxParam : c.minParam();
		new_p = c.evaluate(max_t);
		path->push_back(geoVEC3F(std::array < float, 3 > {{new_p[0], new_p[1], z}}));
		for (int i = c._knot.size()-1; i >= 0; i--){
			if (min_t < c._knot[i] && c._knot[i] < max_t){
				new_p = c.evaluate(c._knot[i]);
				path->push_back(geoVEC3F(std::array < float, 3 > {{new_p[0], new_p[1], z}}));
			}
		}
		new_p = c.evaluate(min_t);
		path->push_back(geoVEC3F(std::array < float, 3 > {{new_p[0], new_p[1], z}}));
	}
}



////cuts along an offset curve between minParam and maxParam at the specified height
//void toolpath_base::cut_curve_const_z(std::vector<geoVEC3F>* path, float tool_diameter, geoCURVE2F& c, float z, double minParam, double maxParam, float resolution)
//{
//	//bspline::offsetCurve2f oc(&c, 0);
//	cut_curve_const_z(path, tool_diameter, c, z, minParam, maxParam, resolution);
//}

//cuts along an offset curve between minParam and maxParam at the specified height.  Path passes through start and finish points
void toolpath_base::cut_curve_const_z(std::vector<geoVEC3F>* path, float tool_diameter, geoCURVE2F& c, float z, double minParam, double maxParam, float resolution)
{
	//c.SetOffset(c.GetOffset() + tool_diameter / 2);

	double min_t, max_t, t,kt;
	geoVEC2F new_p2;
	geoVEC3F new_p3;
	geoVEC2F tangent;

	min_t = (c.minParam() < minParam) ? minParam : c.minParam();
	min_t = (c.maxParam() > minParam) ? minParam : c.maxParam();
	max_t = (c.minParam() < maxParam) ? maxParam : c.minParam();
	max_t = (c.maxParam() > maxParam) ? maxParam : c.maxParam();

	//int prev_span = bspline::find_span(c.GetBase()->lKnot(),c.GetBase()->getKnotVector(),min_t);

	int n = (int)(c.length(min_t, max_t) / resolution);

	for (int i = 0; i <= n; i++){

		double t = min_t * (1 - (double)i / n) + max_t * (double)i / n;
		if (t < c.minParam()){ t = c.minParam(); }
		if (t > c.maxParam()){ t = c.maxParam(); }

		new_p2 = c.evaluate(t);
		new_p3[0] = new_p2[0];
		new_p3[1] = new_p2[1];
		new_p3[2] = z;

		if (i == 0 || i == n){
			//always add points at min_t and max_t so endpoints of line are exact
			cut_curve_add_point(path, new_p3, tangent, 0);
		}
		else{
			//only add intermediate points if deviate from straight line sufficiently
			cut_curve_add_point(path, new_p3, tangent, 0.02);
		}
	}
}


void toolpath_base::cut_curve_add_point(std::vector<geoVEC3F>* path, geoVEC3F new_p, geoVEC2F& tangent, float resolution){
	int size = path->size();

	if (size > 1)
	{
		//float ux = path->operator[](size - 1).x - path->operator[](size - 2).x;
		//float uy = path->operator[](size - 1).y - path->operator[](size - 2).y;
		float vx = new_p[0] - path->operator[](size - 1)[0];
		float vy = new_p[1] - path->operator[](size - 1)[1];

		//only add point if its perpendicular distance to line through previous line segment is greater than 0.1mm
		//always add point if it is first or last in curve
		//float d = abs(ux * vy - uy * vx) / sqrtf(ux * ux + uy * uy) > 0.1;
		float d = abs(tangent[0] * vy - tangent[1] * vx);
		if (d >= resolution){
			tangent[0] = new_p[0] - path->back()[0];
			tangent[1] = new_p[1] - path->back()[1];
			tangent.normalise();
			path->push_back(new_p);
		}
		else{
			path->operator[](size - 1) = new_p;
		}
	}
	else if(size == 1){
		tangent[0] = new_p[0] - path->back()[0];
		tangent[1] = new_p[1] - path->back()[1];
		tangent.normalise();
		path->push_back(new_p);
	}
	else{
		path->push_back(new_p);
	}
}


void toolpath_base::offset_curve(std::vector<geoVEC3F>* path, float tool_diameter,bool deep_corners)
{
	int j,k;
	geoVEC2D param;
	std::vector<geoVEC3F> offset_path;
	offset_path.reserve(path->size());

	for (int i = 0; i < path->size()-2; i++){
		j = (i + 1) % path->size();
		k = (i + 2) % path->size();

		geoVEC3F v1=path->operator[](i);
		geoVEC3F v2=path->operator[](j);
		geoVEC3F v3=path->operator[](k);

		geoVEC2F n1(std::array < float, 2 > {{v1[1] - v2[1], v2[0] - v1[0]}});
		n1.normalise();
		n1 *= tool_diameter * 0.5;
		geoVEC2F n2(std::array < float, 2 > {{v2[1] - v3[1], v3[0] - v2[0]}});
		n2.normalise();
		n2 *= tool_diameter * 0.5;

		geoVEC2F p1(std::array < float, 2 > {{v1[0] + n1[0], v1[1] + n1[1]}});
		geoVEC2F p2(std::array < float, 2 > {{v2[0] + n1[0], v2[1] + n1[1]}});
		geoVEC2F q1(std::array < float, 2 > {{v2[0] + n2[0], v2[1] + n2[1]}});
		geoVEC2F q2(std::array < float, 2 > {{v3[0] + n2[0], v3[1] + n2[1]}});

		bool intersect = geometry::IntersectionParam(p1, p2, q1, q2, param, false);

		if (i == 0){
			offset_path.push_back(geoVEC3F(std::array < float, 3 > {{p1[0], p1[1], v1[2]}}));
		}
		
		geoVEC3F c;
		c[0] = param[0] * p1[0] + (1 - param[0]) * p2[0];
		c[1] = param[0] * p1[1] + (1 - param[0]) * p2[1];
		c[2] = v2[2];
		offset_path.push_back(c);

		//intersect==true if we are inside rather than outside a corner
		if (deep_corners && intersect == true){
			geoVEC2F bisector(std::array < float, 2 > {{c[0] - v2[0], c[1] - v2[1]}});
			bisector.normalise();
			bisector *= tool_diameter * 0.5;

			geoVEC3F v(std::array < float, 3 > {{v2[0] + bisector[0], v2[1] + bisector[1], v2[2]}});
			offset_path.push_back(v);
			offset_path.push_back(c);
		}
				
		if (k == path->size() - 1){
			offset_path.push_back(geoVEC3F(std::array < float, 3 > {{q2[0], q2[1], v3[2]}}));
		}
	}
	path->clear();
	for (int i = 0; i < offset_path.size(); i++){
		path->push_back(offset_path[i]);
	}

}


void toolpath_base::join_curves_external(std::vector<geoVEC3F>& path1, std::vector<geoVEC3F>& path2)
{
	if (path1.size()>1 && path2.size() > 1)
	{
		int l;
		l = path1.size() - 1;
		geoVEC2F p1; p1 << path1[l];
		geoVEC2F p2; p2 << path1[l-1];

		geoVEC2F q1; q1 << path2[0];
		geoVEC2F q2; q2 << path2[1];

		geoVEC2D param;
		bool intersect = geometry::IntersectionParam(p1, p2, q1, q2, param, false);

		geoVEC3F new_p;
		new_p[0] = param[0] * p1[0] + (1 - param[0]) * p2[0];
		new_p[1] = param[0] * p1[1] + (1 - param[0]) * p2[1];
		new_p[2] = path1.back()[2];

		path1.push_back(new_p);
	}
	path1.insert(path1.end(), path2.begin(), path2.end());
	//if (path1.size()>1 && path2.size() > 1)
	//{
	//	int l;
	//	l = path1.size() - 1;
	//	geoVEC2F p1(path1[l].x, path1[l].y);
	//	geoVEC2F p2(path1[l - 1].x, path1[l - 1].y);

	//	geoVEC2F q1(path2[0].x, path2[0].y);
	//	geoVEC2F q2(path2[1].x, path2[1].y);

	//	geoVEC2D param;
	//	bool intersect = bspline::IntersectionParam(p1, p2, q1, q2, &param, false);

	//	path1.push_back(VEC3F(param.x * p1.x + (1 - param.x) * p2.x, param.x * p1.y + (1 - param.x) * p2.y, path1.back().z));
	//}
	//path1.insert(path1.end(), path2.begin(), path2.end());
}

void toolpath_base::close_curve_external(std::vector<geoVEC3F>& path)
{
	if (path.size()>3)
	{
		int l;
		l = path.size() - 1;
		geoVEC2F p1; p1 << path[0];// (path[0].x, path[0].y);
		geoVEC2F p2; p2 << path[1]; // (path[1].x, path[1].y);

		geoVEC2F q1; q1 << path[l];// (path[l].x, path[l].y);
		geoVEC2F q2; q2 << path[l - 1]; // (path[l - 1].x, path[l - 1].y);

		geoVEC2D param;
		bool intersect = geometry::IntersectionParam(p1, p2, q1, q2, param, false);

//		VEC3F new_p(param.x * p1.x + (1 - param.x) * p2.x, param.x * p1.y + (1 - param.x) * p2.y, path.back().z);
		geoVEC3F new_p;
		new_p[0] = param[0] * p1[0] + (1 - param[0]) * p2[0];
		new_p[1] = param[0] * p1[1] + (1 - param[0]) * p2[1];
		new_p[2] = path.back()[2];

		path.push_back(new_p);
		path.push_back(new_p);
		for (int i = path.size() - 2; i > 0; i--){
			path[i] = path[i - 1];
		}
		path[0] = new_p;
	}
}





void toolpath_base::trim_corner_block(geoCURVE2F& c_bout, geoCURVE2F& bout, geoCURVE2F& block, float block_top_z, float block_bottom_z, float step_z, float resolution)
{
	float error = 1E-9;
	std::vector<geoVEC3F> layer;

	float safe_z = parameters["safe_z"];

	geometry::vector<double, 2> param1, param2;
	geometry::bspline::IntersectParam(c_bout, block, error, param1, false);
	IntersectParam(c_bout, bout, error, param2, false);

	param2[0] = (param2[0] > param1[0]) ? c_bout.maxParam() : c_bout.minParam();

	cut_curve_const_z(&layer, tool->diameter, c_bout, 0, param1[0], param2[0], resolution);

	points.push_back(geoVEC3F(std::array < float, 3 > {{ layer[0][0], layer[0][1], safe_z }}));
	float z = block_top_z;
	bool even = true;
	while (z > 	block_bottom_z)
	{
		z -= step_z;
		if (z < block_bottom_z){ z = block_bottom_z; }
		if (even){
			for (int i = 0; i < layer.size(); i++){
				points.push_back(geoVEC3F(std::array < float, 3 > {{layer[i][0], layer[i][1], z}}));
				even = false;
			}
		}
		else{
			for (int i = layer.size() - 1; i >= 0; i--){
				points.push_back(geoVEC3F(std::array < float, 3 > {{layer[i][0], layer[i][1], z}}));
				even = true;
			}
		}
	}
	int n = layer.size() - 1;
	points.push_back(geoVEC3F(std::array < float, 3 > {{points.back()[0], points.back()[1], safe_z}}));
}

void toolpath_base::trim_end_block(geoCURVE2F& rib, geoCURVE2F& centre_line, geoCURVE2F& block, float block_top_z, float block_bottom_z, float step_z, float resolution)
{
	float error = 1E-9;
	std::vector<geoVEC3F> layer;
	float safe_z = parameters["safe_z"];

	geoVEC2D param1, param2, param3;

	IntersectParam(block, centre_line, error, param1, false);

	geoCURVE2F block1 = block; block1.trim_curve(block.minParam(), param1[0]);
	geoCURVE2F block2 = block; block2.trim_curve(param1[0], block.maxParam());

	IntersectParam(rib, block1, error, param2, false);
	IntersectParam(rib, block2, error, param3, false);

	cut_curve_const_z(&layer, tool->diameter + 1.0, rib, 0, param2[0], param3[0], resolution);

	points.push_back(geoVEC3F(std::array < float, 3 > {{layer[0][0], layer[0][1], safe_z}}));
	float z = block_top_z;
	bool even = true;
	while (z > 	block_bottom_z)
	{
		z -= step_z;
		if (z < block_bottom_z){ z = block_bottom_z; }

		if (even){
			for (int i = 0; i < layer.size(); i++){
				points.push_back(geoVEC3F(std::array < float, 3 > {{layer[i][0], layer[i][1], z}}));
			}
			even = false;
		}
		else{
			for (int i = layer.size() - 1; i >= 0; i--){
				points.push_back(geoVEC3F(std::array < float, 3 > {{layer[i][0], layer[i][1], z}}));
			}
			even = true;
		}
	}

	//shave in to final edge
	double u, v;
	if (even) { u = param2[0]; v = param3[0]; }
	else{ u = param3[0]; v = param2[0]; }

	cut_curve_const_z(&points, tool->diameter + 0.5, rib, block_bottom_z, u, v, resolution);
	cut_curve_const_z(&points, tool->diameter + 0.0, rib, block_bottom_z, v, u, resolution);

	points.push_back(geoVEC3F(std::array < float, 3 > {{points.back()[0], points.back()[1], safe_z}}));

}

/*
This method calculates tool path over surface.  It takes passes parallel to x axis, with each pass incrementing y coordinate.

Each pass parallel to x axis is divided into steps of "step_x" mm. At each (x,y) coordinate, the z coordinate is calculated 
to be the minimum z value such that the tool tip does not intersect any triangles of the mesh.

The method assumes a ball-nose tool is used.

CAUTION: step sizes must be small enough so path does not cut significantly into a convex surface
*/
void toolpath_base::finish_surface_scanning_stl(std::vector<geoVEC3F>* path, float tool_diameter, float step_x, float step_y, float minimum_z, geometry::mesh3f& mesh)
{
	if (step_y > tool_diameter * 0.5) return;
	
	mesh.build_nodes(50, 50);





	float min_x = FLT_MAX, max_x = FLT_MIN, min_y = FLT_MAX, max_y = FLT_MIN, min_z = FLT_MAX, max_z = FLT_MIN;
	
	//1. find ranges of stl
	min_y = mesh.node_y.front() - 0.5 * tool_diameter;
	max_y = mesh.node_y.back() + 0.5 * tool_diameter;


	//2. move over stl taking steps in y direction
	int steps_y = (int)((max_y - min_y) / step_y + 1);

	float x, x_plus, x_minus, y, y_plus, y_minus, z;

	std::vector<int> node_cols;

	for (int j = 0; j <= steps_y; j++)
	{
		
		//2.1 add all triangles that could overlap current y strip to strip vector
		y = (1 - (float) j / steps_y) * min_y + (float)j / steps_y * max_y;
		y_plus = y + 0.5 * tool_diameter;
		y_minus = y - 0.5 * tool_diameter;

		node_cols.clear();

		for (int col = 0; col < mesh.node_y.size()-1; col++)
		{
			if (y_plus < mesh.node_y[col] || y_minus > mesh.node_y[col + 1]) continue;
			
			node_cols.push_back(col); 
			min_x = FLT_MAX; 
			max_x = FLT_MIN;
			for (int row = 0; row < mesh.node_x.size() - 1; row++)
			{
				if (mesh.nodes.GetPoint(row, col).size() > 0){
					min_x = fminf(min_x, mesh.node_x[row]);
					max_x = fmaxf(max_x, mesh.node_x[row+1]);
				}
			}
			
		}
		min_x -= 0.5 * tool_diameter;
		max_x += 0.5 * tool_diameter;

		//2.2 move over strip in x direction 
		int steps_x = (int)((max_x - min_x) / step_x + 1);

		geoVEC2F p1, p2, p3;

		for (int i = 0; i < steps_x; i++)
		{
			x = (1 - (float)i / steps_x) * min_x + (float)i / steps_x * max_x;
			x_plus = x + 0.5 * tool_diameter;
			x_minus = x - 0.5 * tool_diameter;

			z = minimum_z;

			for (int row = 0; row < mesh.node_x.size() - 1; row++)
			{
				if (x_plus < mesh.node_x[row] || x_minus > mesh.node_x[row + 1]) continue;
				
				for (int col : node_cols)
				{
					for (int e : mesh.nodes.GetPoint(row, col))
					{
						geoVEC2F p1; p1 << mesh.get_vertex(e, 0);
						geoVEC2F p2; p2 << mesh.get_vertex(e, 1);
						geoVEC2F p3; p3 << mesh.get_vertex(e, 2);

						if (geometry::interior_triangle(geoVEC2F(std::array < float, 2 > {{x, y}}), (float) 0.5 * tool_diameter, p1, p2, p3))
						{
							//z = fmaxf(z, mesh.get_vertex(e, 0)[2]);
							//z = fmaxf(z, mesh.get_vertex(e, 1)[2]);
							//z = fmaxf(z, mesh.get_vertex(e, 2)[2]);

							geoVEC3F q1 = mesh.get_vertex(e, 0);
							geoVEC3F q2 = mesh.get_vertex(e, 1);
							geoVEC3F q3 = mesh.get_vertex(e, 2);

							float r = (float) 0.5 * tool_diameter;

							//find heights at which tool tip touches vertices of triangle.  If no touches, DBL_MIN is returned
							z = fmaxf(z, geometry::min_height_sphere_on_point(x, y, r, q1));
							z = fmaxf(z, geometry::min_height_sphere_on_point(x, y, r, q2));
							z = fmaxf(z, geometry::min_height_sphere_on_point(x, y, r, q3));

							//find heights at which tool tip touches edges of triangle.  If no touches, DBL_MIN is returned
							z = fmaxf(z, geometry::min_height_sphere_on_line(x, y, r, q1, q2));
							z = fmaxf(z, geometry::min_height_sphere_on_line(x, y, r, q2, q3));
							z = fmaxf(z, geometry::min_height_sphere_on_line(x, y, r, q3, q1));

							//find height at which tool tip touches triangle
							z = fmaxf(z, geometry::min_height_sphere_on_triangle(x, y, r, q1, q2, q3));
						}

					}
				}
			}
			if (z > DBL_MIN)
			{
				path->push_back(geoVEC3F(std::array < float, 3 > {{x, y, z}}));
			}
		}


	}


}



/*
Method calculates finishing toolpath by tracing over mesh in x,y plane.  For each row, many points are calculated which lie on each 
triangle in regin of tool tip.  The z coordinate is calculated to be the minimum z such that the tool is above all points.

Enough points must be created so that the tool tip doesn't significanty dip into the surface if between points
*/
void toolpath_base::finish_surface_point_cloud(std::vector<geoVEC3F>* path, float tool_diameter, float step_x, float step_y, float minimum_z, geometry::mesh3f& mesh)
{

	int n = 20;
	float u;
	float v;

	std::vector<geoVEC3F> dense_points;
	for (int e = 0; e < mesh.elements.size(); e++)
	{
		geoVEC3F p1 = mesh.get_vertex(e, 0);
		geoVEC3F p2 = mesh.get_vertex(e, 1);
		geoVEC3F p3 = mesh.get_vertex(e, 2);
		geoVEC3F new_p;
		for (int i = 0; i <= n; i++)
		{
			u = (float)i / n;
			for (int j = 0; j <= n - i; j++)
			{
				v = (double)j / n;
				new_p = (1 - u - v) * p1 + u * p2 + v * p3;
				dense_points.push_back(new_p);
			}
		}

	}
	std::sort(dense_points.begin(), dense_points.begin(), [](geoVEC3F const& a, geoVEC3F const& b){return a[1]<b[1]; });


}


void toolpath_base::rough_surface_grid(std::vector<geoVEC3F>* path, float tool_diameter, float step_x, float step_y, float minimum_z, float safe_z, float margin_z, geometry::mesh3f& mesh)
{
	
	if (step_y > tool_diameter * 0.5) return;
	if (margin_z < 0) return;

	geoVEC3F min, max;
	mesh.get_mesh_range(min, max);


	//1. construct grid of max heights over mesh
	std::vector<float> grid_x;
	std::vector<float> grid_y;

	float x = min[0];
	do
	{
		grid_x.push_back(x);
		x += step_x;
	} while (x < max[0] + step_x);

	float y = min[1];
	do
	{
		grid_y.push_back(y);
		y += step_y;
	} while (y < max[1] + step_y);

	geometry::lattice<float> grid_mesh_heights(grid_x.size() - 1, grid_y.size() - 1);

	for (int g = 0; g < grid_mesh_heights.data.size(); g++)
	{
		grid_mesh_heights.data[g] = FLT_MIN;
	}

	for (int e = 0; e < mesh.elements.size(); e++){

		mesh.get_element_range(e, min, max);
		for (int i = 0; i < grid_x.size() - 1; i++){
			if (min[0] <= grid_x[i + 1] && max[0] >= grid_x[i])
			{
				for (int j = 0; j < grid_y.size() - 1; j++){
					if (min[1] <= grid_y[j + 1] && max[1] >= grid_y[j])
					{
						grid_mesh_heights.SetPoint(i, j, max[2] + margin_z);
					}
				}
			}
		}
	}


	//2. step over grid using max heights to set tool height
	float y_min, y_max, x_min, x_max, z;
	int start, end, incr;
	int row_width, col_width;

	row_width = (int)(tool_diameter * 0.5) / step_x + 1;
	col_width = (int)(tool_diameter * 0.5) / step_y + 1;

	bool at_safe_z = true;
	//bool valid_row;

	for (int j = 0; j < grid_y.size(); j++)
	{
		y = grid_y[j];

		start = (j % 2 == 0) ? 0 : grid_x.size() - 2;
		end = (j % 2 == 0) ? grid_x.size() - 1 : -1;
		incr = (j % 2 == 0) ? 1 : -1;

		for (int i = start; i != end; i+=incr)
		{
			
			x = grid_x[i];

			z = FLT_MIN;
			for (int j0 = fmaxf(0, j - col_width); j0 < fminf(grid_mesh_heights.cols, j + col_width); j0++){
				for (int i0 = fmaxf(0, i - row_width - 1); i0 < fminf(grid_mesh_heights.rows, i + row_width + 1); i0++){
					z = fmaxf(z, grid_mesh_heights.GetPoint(i0, j0));
				}
			}

			if (z > FLT_MIN){
				if (at_safe_z == true){
					path->push_back(geoVEC3F(std::array < float, 3 > {{x, y, safe_z}}));
					at_safe_z = false;
				}
				path->push_back(geoVEC3F(std::array < float, 3 > {{x, y, z}}));
			}
		}
	}
	geoVEC3F p = path->back();
	path->push_back(geoVEC3F(std::array < float, 3 > {{p[0], p[1], safe_z}}));


}




