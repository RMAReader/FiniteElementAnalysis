
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

		out << "G01";
		out << "X" << fixed << setprecision(3) << points[0].x;
		out << "Y" << fixed << setprecision(3) << points[0].y;
		out << "F" << fixed << setprecision(0) << parameters["max_feedrate"];
		out << endl;

		for (int i = 0; i < points.size(); i++)
		{
			out << "G01";
			out << "X" << fixed << setprecision(3) << points[i].x;
			out << "Y" << fixed << setprecision(3) << points[i].y;
			out << "Z" << fixed << setprecision(3) << points[i].z;
			out << "F" << fixed << setprecision(0) << parameters["xy_feedrate"];
			out << endl;
		}

		out << "M05"; //spindle off
		out << "M30"; //end program
		out << endl;

		out.close();
	}

}


void toolpath_base::bore_hole(circle2f& c, float safe_z, float max_z, float min_z, float step_z)
{
	if (tool == nullptr){ exit; }
	if (tool->diameter < 0){ exit; }


	float radius = (c.diameter - tool->diameter) / 2;
	if (radius > 0)
	{

		int n = 100;
		int j = 0;
		double PI = 3.141592653589793238463;
		float z = max_z;
		float theta = 0;

		//move to safe z and centre of circle
		points.push_back(point3f(c.centre.x, c.centre.y, safe_z));
		points.push_back(point3f(c.centre.x, c.centre.y, z));

		//spiral down to required depth
		while (z >= min_z)
		{
			theta = PI * 2 * j / n;
			z -= step_z / n;
			if (z >= min_z)
			{
				point3f p;
				p.x = c.centre.x + radius * sinf(theta);
				p.y = c.centre.y + radius * cosf(theta);
				p.z = z;
				points.push_back(p);
				j++;
			}
		}

		//do one full circle at required depth
		for (int k = j; k <= j + n; k++)
		{
			theta = PI * 2 * k / n;
			point3f p;
			p.x = c.centre.x + radius * sinf(theta);
			p.y = c.centre.y + radius * cosf(theta);
			p.z = min_z;
			points.push_back(p);
		}

		//move to centre of circle and safe z
		points.push_back(point3f(c.centre.x, c.centre.y, min_z));
		points.push_back(point3f(c.centre.x, c.centre.y, safe_z));
	}
}



void toolpath_base::cut_polyline_const_z(std::vector<VEC3F>* path, bspline::curve<bspline::vec2<float>>& c, float z, double minParam, double maxParam)
{
	double min_t, max_t;
	VEC2F new_p;

	if (minParam < maxParam){
		min_t = (c.minParam() < minParam) ? minParam : c.minParam();
		max_t = (c.maxParam() > maxParam) ? maxParam : c.maxParam();
		new_p = c.evaluate(min_t);
		path->push_back(VEC3F(new_p.x,new_p.y,z));
		for (int i = 0; i < c.lKnot(); i++){
			if (min_t < c.getKnot(i) && c.getKnot(i) < max_t){
				new_p = c.evaluate(c.getKnot(i));
				path->push_back(VEC3F(new_p.x, new_p.y, z));
			}
		}
		new_p = c.evaluate(max_t);
		path->push_back(VEC3F(new_p.x, new_p.y, z));
	}
	else{
		max_t = (c.maxParam() > minParam) ? minParam : c.maxParam();
		min_t = (c.minParam() < maxParam) ? maxParam : c.minParam();
		new_p = c.evaluate(max_t);
		path->push_back(VEC3F(new_p.x, new_p.y, z));
		for (int i = c.lKnot()-1; i >= 0; i--){
			if (min_t < c.getKnot(i) && c.getKnot(i) < max_t){
				new_p = c.evaluate(c.getKnot(i));
				path->push_back(VEC3F(new_p.x, new_p.y, z));
			}
		}
		new_p = c.evaluate(min_t);
		path->push_back(VEC3F(new_p.x, new_p.y, z));
	}
}



//cuts along an offset curve between minParam and maxParam at the specified height
void toolpath_base::cut_curve_const_z(std::vector<VEC3F>* path, float tool_diameter, CURVE2F& c, float z, double minParam, double maxParam, float resolution)
{
	bspline::offsetCurve2f oc(&c, 0);
	cut_curve_const_z(path, tool_diameter, oc, z, minParam, maxParam, resolution);
}

//cuts along an offset curve between minParam and maxParam at the specified height.  Path passes through start and finish points
void toolpath_base::cut_curve_const_z(std::vector<VEC3F>* path, float tool_diameter, bspline::offsetCurve2f& c, float z, double minParam, double maxParam, float resolution)
{
	c.SetOffset(c.GetOffset() + tool_diameter / 2);

	double min_t, max_t, t,kt;
	VEC2F new_p;
	VEC2F tangent;

	min_t = (c.GetBase()->minParam() < minParam) ? minParam : c.GetBase()->minParam();
	min_t = (c.GetBase()->maxParam() > minParam) ? minParam : c.GetBase()->maxParam();
	max_t = (c.GetBase()->minParam() < maxParam) ? maxParam : c.GetBase()->minParam();
	max_t = (c.GetBase()->maxParam() > maxParam) ? maxParam : c.GetBase()->maxParam();

	//int prev_span = bspline::find_span(c.GetBase()->lKnot(),c.GetBase()->getKnotVector(),min_t);

	int n = (int)(c.GetBase()->length(min_t, max_t) / resolution);

	for (int i = 0; i <= n; i++){

		double t = min_t * (1 - (double)i / n) + max_t * (double)i / n;
		if (t < c.GetBase()->minParam()){ t = c.GetBase()->minParam(); }
		if (t > c.GetBase()->maxParam()){ t = c.GetBase()->maxParam(); }

		new_p= c.evaluate(t);
		if (i == 0 || i == n){
			//always add points at min_t and max_t so endpoints of line are exact
			cut_curve_add_point(path, VEC3F(new_p.x, new_p.y, z),tangent, 0);
		}
		else{
			//only add intermediate points if deviate from straight line sufficiently
			cut_curve_add_point(path, VEC3F(new_p.x, new_p.y, z),tangent, 0.02);
		}
	}
}


void toolpath_base::cut_curve_add_point(std::vector<VEC3F>* path, VEC3F new_p, VEC2F& tangent, float resolution){
	int size = path->size();

	if (size > 1)
	{
		//float ux = path->operator[](size - 1).x - path->operator[](size - 2).x;
		//float uy = path->operator[](size - 1).y - path->operator[](size - 2).y;
		float vx = new_p.x - path->operator[](size - 1).x;
		float vy = new_p.y - path->operator[](size - 1).y;

		//only add point if its perpendicular distance to line through previous line segment is greater than 0.1mm
		//always add point if it is first or last in curve
		//float d = abs(ux * vy - uy * vx) / sqrtf(ux * ux + uy * uy) > 0.1;
		float d = abs(tangent.x * vy - tangent.y * vx);
		if (d >= resolution){
			tangent.x = new_p.x - path->back().x;
			tangent.y = new_p.y - path->back().y;
			tangent *= 1 / sqrtf(tangent.x * tangent.x + tangent.y * tangent.y);
			path->push_back(new_p);
		}
		else{
			path->operator[](size - 1) = new_p;
		}
	}
	else if(size == 1){
		tangent.x = new_p.x - path->back().x;
		tangent.y = new_p.y - path->back().y;
		tangent *= 1 / sqrtf(tangent.x * tangent.x + tangent.y * tangent.y);
		path->push_back(new_p);
	}
	else{
		path->push_back(new_p);
	}
}


void toolpath_base::offset_curve(std::vector<VEC3F>* path, float tool_diameter,bool deep_corners)
{
	int j,k;
	VEC2D param;
	std::vector<VEC3F> offset_path;
	offset_path.reserve(path->size());

	for (int i = 0; i < path->size()-2; i++){
		j = (i + 1) % path->size();
		k = (i + 2) % path->size();

		VEC3F v1=path->operator[](i);
		VEC3F v2=path->operator[](j);
		VEC3F v3=path->operator[](k);

		VEC2F n1(v1.y - v2.y, v2.x - v1.x);
		n1 *= tool_diameter * 0.5 / sqrtf(n1.x * n1.x + n1.y * n1.y);
		VEC2F n2(v2.y - v3.y, v3.x - v2.x);
		n2 *= tool_diameter * 0.5 / sqrtf(n2.x * n2.x + n2.y * n2.y);

		VEC2F p1(v1.x + n1.x, v1.y + n1.y);
		VEC2F p2(v2.x + n1.x, v2.y + n1.y);
		VEC2F q1(v2.x + n2.x, v2.y + n2.y);
		VEC2F q2(v3.x + n2.x, v3.y + n2.y);

		bool intersect = bspline::IntersectionParam(p1, p2, q1, q2, &param, false);

		if (i == 0){
			offset_path.push_back(VEC3F(p1.x,p1.y,v1.z));
		}
		
		VEC3F c(param.x * p1.x + (1 - param.x) * p2.x, param.x * p1.y + (1 - param.x) * p2.y, v2.z);
		offset_path.push_back(c);

		//intersect==true if we are inside rather than outside a corner
		if (deep_corners && intersect == true){
			VEC2F bisector(c.x - v2.x, c.y - v2.y);
			bisector *= 1 / sqrtf(bisector.x*bisector.x + bisector.y*bisector.y);
			
			VEC3F v(v2.x + bisector.x * tool_diameter * 0.5, v2.y + bisector.y * tool_diameter * 0.5, v2.z);
			offset_path.push_back(v);
			offset_path.push_back(c);
		}
				
		if (k == path->size() - 1){
			offset_path.push_back(VEC3F(q2.x,q2.y,v3.z));
		}
	}
	path->clear();
	for (int i = 0; i < offset_path.size(); i++){
		path->push_back(offset_path[i]);
	}

}


void toolpath_base::join_curves_external(std::vector<VEC3F>& path1, std::vector<VEC3F>& path2)
{
	if (path1.size()>1 && path2.size() > 1)
	{
		int l;
		l = path1.size() - 1;
		VEC2F p1(path1[l].x, path1[l].y);
		VEC2F p2(path1[l - 1].x, path1[l - 1].y);

		VEC2F q1(path2[0].x, path2[0].y);
		VEC2F q2(path2[1].x, path2[1].y);

		VEC2D param;
		bool intersect = bspline::IntersectionParam(p1, p2, q1, q2, &param, false);

		path1.push_back(VEC3F(param.x * p1.x + (1 - param.x) * p2.x, param.x * p1.y + (1 - param.x) * p2.y, path1.back().z));
	}
	path1.insert(path1.end(), path2.begin(), path2.end());
}

void toolpath_base::close_curve_external(std::vector<VEC3F>& path)
{
	if (path.size()>3)
	{
		int l;
		l = path.size() - 1;
		VEC2F p1(path[0].x, path[0].y);
		VEC2F p2(path[1].x, path[1].y);

		VEC2F q1(path[l].x, path[l].y);
		VEC2F q2(path[l-1].x, path[l-1].y);

		VEC2D param;
		bool intersect = bspline::IntersectionParam(p1, p2, q1, q2, &param, false);

		VEC3F new_p(param.x * p1.x + (1 - param.x) * p2.x, param.x * p1.y + (1 - param.x) * p2.y, path.back().z);
		
		path.push_back(new_p);
		path.push_back(new_p);
		for (int i = path.size() - 2; i > 0; i--){
			path[i] = path[i - 1];
		}
		path[0] = new_p;
	}
}





void toolpath_base::trim_corner_block(CURVE2F* c_bout, CURVE2F* bout, CURVE2F* block, float block_top_z, float block_bottom_z, float step_z, float resolution)
{
	float error = 1E-9;
	std::vector<VEC3F> layer;

	float safe_z = parameters["safe_z"];

	VEC2D param1, param2;
	IntersectParam(*c_bout, *block, error, &param1, false);
	IntersectParam(*c_bout, *bout, error, &param2, false);

	param2.x = (param2.x > param1.x) ? c_bout->maxParam() : c_bout->minParam();

	cut_curve_const_z(&layer, tool->diameter, *c_bout, 0, param1.x, param2.x, resolution);

	points.push_back(VEC3F(layer[0].x, layer[0].y, safe_z));
	float z = block_top_z;
	bool even = true;
	while (z > 	block_bottom_z)
	{
		z -= step_z;
		if (z < block_bottom_z){ z = block_bottom_z; }
		if (even){
			for (int i = 0; i < layer.size(); i++){
				points.push_back(VEC3F(layer[i].x, layer[i].y, z));
				even = false;
			}
		}
		else{
			for (int i = layer.size() - 1; i >= 0; i--){
				points.push_back(VEC3F(layer[i].x, layer[i].y, z));
				even = true;
			}
		}
	}
	int n = layer.size() - 1;
	points.push_back(VEC3F(points.back().x, points.back().y, safe_z));
}

void toolpath_base::trim_end_block(CURVE2F* rib, CURVE2F* centre_line, CURVE2F* block, float block_top_z, float block_bottom_z, float step_z, float resolution)
{
	float error = 1E-9;
	std::vector<VEC3F> layer;
	float safe_z = parameters["safe_z"];

	VEC2D param1, param2, param3;

	IntersectParam(*block, *centre_line, error, &param1, false);

	CURVE2F* block1 = block->trim_curve(block->minParam(), param1.x);
	CURVE2F* block2 = block->trim_curve(param1.x, block->maxParam());

	IntersectParam(*rib, *block1, error, &param2, false);
	IntersectParam(*rib, *block2, error, &param3, false);

	cut_curve_const_z(&layer, tool->diameter + 1.0, *rib, 0, param2.x, param3.x, resolution);

	points.push_back(VEC3F(layer[0].x, layer[0].y, safe_z));
	float z = block_top_z;
	bool even = true;
	while (z > 	block_bottom_z)
	{
		z -= step_z;
		if (z < block_bottom_z){ z = block_bottom_z; }

		if (even){
			for (int i = 0; i < layer.size(); i++){
				points.push_back(VEC3F(layer[i].x, layer[i].y, z));
			}
			even = false;
		}
		else{
			for (int i = layer.size() - 1; i >= 0; i--){
				points.push_back(VEC3F(layer[i].x, layer[i].y, z));
			}
			even = true;
		}
	}

	//shave in to final edge
	double u, v;
	if (even) { u = param2.x; v = param3.x; }
	else{ u = param3.x; v = param2.x; }

	cut_curve_const_z(&points, tool->diameter + 0.5, *rib, block_bottom_z, u, v, resolution);
	cut_curve_const_z(&points, tool->diameter + 0.0, *rib, block_bottom_z, v, u, resolution);

	points.push_back(VEC3F(points.back().x, points.back().y, safe_z));

}


void toolpath_base::rough_surface_scanning_stl(std::vector<VEC3F>* path, float tool_diameter, float step_x, float step_y, float minimum_z, std::vector<TRIANGLE3F>* triangles)
{
	if (step_y > tool_diameter * 0.5) return;
	
	float min_x = FLT_MAX, max_x = FLT_MIN, min_y = FLT_MAX, max_y = FLT_MIN, min_z = FLT_MAX, max_z = FLT_MIN;
	
	//1. find ranges of stl
	for (auto t = triangles->begin(); t != triangles->end(); ++t){
		min_y = fminf(min_y, t->min_y());
		max_y = fmaxf(max_y, t->max_y());
	}
	min_y -= 0.5 * tool_diameter;
	max_y += 0.5 * tool_diameter;


	//2. move over stl taking steps in y direction
	int steps_y = (int)((max_y - min_y) / step_y + 1);

	float x, x_plus, x_minus, y, y_plus, y_minus, z;
	std::vector<TRIANGLE3F*> strip;

	for (int j = 0; j <= steps_y; j++)
	{
		
		//2.1 add all triangles that could overlap current y strip to strip vector
		y = (1 - (float) j / steps_y) * min_y + (float)j / steps_y * max_y;
		y_plus = y + 0.5 * tool_diameter;
		y_minus = y - 0.5 * tool_diameter;

		strip.clear();
		for (auto t = triangles->begin(); t != triangles->end(); ++t)
		{
			if (t->max_y() > y_minus && t->min_y() < y_plus)
			{  
				strip.push_back(t._Ptr);
				min_x = fminf(min_x, t->min_x());
				max_x = fmaxf(max_x, t->max_x());
			}
		}
		min_x -= 0.5 * tool_diameter;
		max_x += 0.5 * tool_diameter;

		//2.2 move over strip in x direction 
		int steps_x = (int)((max_x - min_x) / step_x + 1);

		for (int i = 0; i < steps_x; i++)
		{
			x = (1 - (float)i / steps_x) * min_x + (float)i / steps_x * max_x;
			x_plus = x + 0.5 * tool_diameter;
			x_minus = x - 0.5 * tool_diameter;

			z = minimum_z;

			for (auto t = strip.begin(); t != strip.end(); ++t)
			{
				if ((*t)->max_x() > x_minus && (*t)->min_x() < x_plus)
				{
					if (bspline::interior_triangle(VEC2F(x, y), 0.5 * tool_diameter, VEC2F((*t)->p[0].x, (*t)->p[0].y), VEC2F((*t)->p[1].x, (*t)->p[1].y), VEC2F((*t)->p[2].x, (*t)->p[2].y)))
					{
						z = fmaxf(z, (*t)->max_z());
					}
				}
			}
			path->push_back(VEC3F(x, y, z));
		}


	}


}


void toolpath_RibMould::calculate()
{
	//zero z is set to top of rib mould base

	using namespace bspline;

	float safe_z = parameters["safe_z"];
	float mould_top_z = parameters["mould_top_z"];
	float mould_bottom_z = parameters["mould_bottom_z"];
	float step_z = parameters["step_z"];
	

	//bore location holes
	for (int i = 0; i < violin->ribs->rib_mould_locator_holes.size(); i++)
	{
		bore_hole(violin->ribs->rib_mould_locator_holes[i], safe_z, mould_top_z, mould_bottom_z,step_z);
	}

	float resolution = 1.0;


	//cut out central clamping hole
	std::vector<VEC3F> clamp_hole_layer;
	std::unordered_map<std::string, CURVE2F*>& c = violin->ribs->curves;

	cut_curve_const_z(&clamp_hole_layer, tool->diameter, *c["clamp_hole"], 0, c["clamp_hole"]->minParam(), c["clamp_hole"]->maxParam(), resolution);

	points.push_back(VEC3F(clamp_hole_layer[0].x, clamp_hole_layer[0].y, safe_z));
	float z = mould_top_z;
	while (z > mould_bottom_z)
	{
		z -= step_z;
		if (z < mould_bottom_z){ z = mould_bottom_z; }

		for (int i = 0; i < clamp_hole_layer.size(); i++){
			points.push_back(VEC3F(clamp_hole_layer[i].x, clamp_hole_layer[i].y, z));
		}
	}
	int n = clamp_hole_layer.size() - 1;
	points.push_back(VEC3F(clamp_hole_layer[n].x, clamp_hole_layer[n].y, safe_z));



	//cut outline - clockwise from lower bout treble corner
	float error = 1E-4;

	//lower bout
	bool intersect1, intersect2, intersect3, intersect4, intersect5, intersect6;
	bool intersect7,intersect8, intersect9, intersect10, intersect11, intersect12;
	vec2<double> param1, param2, param3, param4;

	intersect1 = IntersectParam(*c["rib_internal_lower_bout"], *c["block_trbl_lower_bout"], error, &param1, false);
	intersect2 = IntersectParam(*c["rib_internal_lower_bout"], *c["centre_line"], error, &param2, false);
	intersect4 = IntersectParam(*c["rib_internal_lower_bout"], *c["block_bass_lower_bout"], error, &param4, false);
	if (intersect1 && intersect2 && intersect4){
		CURVE2F* rib_lower_trbl = c["rib_internal_lower_bout"]->trim_curve(param1.x, param2.x);
		CURVE2F* rib_lower_bass = c["rib_internal_lower_bout"]->trim_curve(param2.x, param4.x);

		intersect2 = IntersectParam(*rib_lower_trbl, *c["block_bottom"], error, &param2, false);
		intersect3 = IntersectParam(*rib_lower_bass, *c["block_bottom"], error, &param3, false);
	}
	//centre bout bass
	vec2<double> param5, param6;
	intersect5 = IntersectParam(*c["rib_internal_centre_bout_bass"], *c["block_bass_lower_bout"], error, &param5, false);
	intersect6 = IntersectParam(*c["rib_internal_centre_bout_bass"], *c["block_bass_upper_bout"], error, &param6, false);

	//upper bout
	vec2<double> param7, param8, param9, param10;

	intersect7 = IntersectParam(*c["rib_internal_upper_bout"], *c["block_bass_upper_bout"], error, &param7, false);
	intersect8 = IntersectParam(*c["rib_internal_upper_bout"], *c["centre_line"], error, &param8, false);
	intersect10 = IntersectParam(*c["rib_internal_upper_bout"], *c["block_trbl_upper_bout"], error, &param10, false);
	if (intersect7 && intersect8 && intersect10){
		CURVE2F* rib_upper_bass = c["rib_internal_upper_bout"]->trim_curve(param7.x, param8.x);
		CURVE2F* rib_upper_trbl = c["rib_internal_upper_bout"]->trim_curve(param8.x, param10.x);

		intersect8 = IntersectParam(*rib_upper_bass, *c["block_neck"], error, &param8, false);
		intersect9 = IntersectParam(*rib_upper_trbl, *c["block_neck"], error, &param9, false);
	}

	//centre bout treble
	vec2<double> param11, param12;
	intersect11 = IntersectParam(*c["rib_internal_centre_bout_trbl"], *c["block_trbl_upper_bout"], error, &param11, false);
	intersect12 = IntersectParam(*c["rib_internal_centre_bout_trbl"], *c["block_trbl_lower_bout"], error, &param12, false);

	if (intersect1 && intersect2 && intersect3 && intersect4 && intersect5 && intersect6 && 
		intersect7 && intersect8 && intersect9 && intersect10 && intersect11 && intersect12)
	{
		std::vector<VEC3F> layer1,layer2,layer3,layer4,layer5,layer6,layer7,layer8,layer9,layer10,layer11,layer12;

		cut_curve_const_z(&layer1, 0, *c["rib_internal_lower_bout"], 0, param1.x, param2.x, resolution);
		cut_polyline_const_z(&layer2, *c["block_bottom"], 0, param2.y, param3.y);
		cut_curve_const_z(&layer3, 0, *c["rib_internal_lower_bout"], 0, param3.x, param4.x, resolution);

		cut_polyline_const_z(&layer4, *c["block_bass_lower_bout"], 0, param4.y, param5.y);
		cut_curve_const_z(&layer5, 0,*c["rib_internal_centre_bout_bass"], 0, param5.x, param6.x,resolution);
		cut_polyline_const_z(&layer6, *c["block_bass_upper_bout"], 0, param6.y, param7.y);

		cut_curve_const_z(&layer7, 0, *c["rib_internal_upper_bout"], 0, param7.x, param8.x, resolution);
		cut_polyline_const_z(&layer8, *c["block_neck"], 0, param8.y, param9.y);
		cut_curve_const_z(&layer9, 0, *c["rib_internal_upper_bout"], 0, param9.x, param10.x, resolution);

		cut_polyline_const_z(&layer10,  *c["block_trbl_upper_bout"], 0, param10.y, param11.y);
		cut_curve_const_z(&layer11, 0, *c["rib_internal_centre_bout_trbl"], 0, param11.x, param12.x, resolution);
		cut_polyline_const_z(&layer12, *c["block_trbl_lower_bout"], 0, param12.y, param1.y);

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

		std::vector<VEC3F> layer;
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

		points.push_back(VEC3F(layer.front().x, layer.front().y, safe_z));
		z = mould_top_z;
		while (z > mould_bottom_z)
		{
			z -= step_z;
			if (z < mould_bottom_z){ z = mould_bottom_z; }

			for (int i = 0; i < layer.size(); i++){
				points.push_back(VEC3F(layer[i].x, layer[i].y, z));
			}
		}
		points.push_back(VEC3F(layer.back().x, layer.back().y, safe_z));
	}
}

void toolpath_TrimBlocksCentreBout::calculate()
{
	using namespace bspline;

	float safe_z = parameters["safe_z"];
	float block_top_z = parameters["block_top_z"];
	float block_bottom_z = parameters["block_bottom_z"];
	float step_z = parameters["step_z"];
	float resolution = 0.1;

	std::unordered_map<std::string, CURVE2F*>& c = violin->ribs->curves;
	points.push_back(VEC3F(0, 0, safe_z));

	trim_corner_block(c["rib_internal_centre_bout_trbl"], c["rib_internal_lower_bout"], c["block_trbl_lower_bout"], block_top_z, block_bottom_z, step_z, resolution);
	trim_corner_block(c["rib_internal_centre_bout_trbl"], c["rib_internal_upper_bout"], c["block_trbl_upper_bout"], block_top_z, block_bottom_z, step_z, resolution);
	trim_corner_block(c["rib_internal_centre_bout_bass"], c["rib_internal_upper_bout"], c["block_bass_upper_bout"], block_top_z, block_bottom_z, step_z, resolution);
	trim_corner_block(c["rib_internal_centre_bout_bass"], c["rib_internal_lower_bout"], c["block_bass_lower_bout"], block_top_z, block_bottom_z, step_z, resolution);

}

void toolpath_TrimBlocksEndBouts::calculate()
{
	using namespace bspline;

	float safe_z = parameters["safe_z"];
	float block_top_z = parameters["block_top_z"];
	float block_bottom_z = parameters["block_bottom_z"];
	float step_z = parameters["step_z"];
	float resolution = 0.1;

	std::unordered_map<std::string, CURVE2F*>& c = violin->ribs->curves;

	points.push_back(VEC3F(0, 0, safe_z));

	trim_corner_block(c["rib_internal_lower_bout"], c["rib_internal_centre_bout_trbl"], c["block_trbl_lower_bout"], block_top_z, block_bottom_z, step_z, resolution);
	trim_end_block(c["rib_internal_lower_bout"], c["centre_line"], c["block_bottom"], 31, 0, 5, 0.1);
	trim_corner_block(c["rib_internal_lower_bout"], c["rib_internal_centre_bout_bass"], c["block_bass_lower_bout"], block_top_z, block_bottom_z, step_z, resolution);

	trim_corner_block(c["rib_internal_upper_bout"], c["rib_internal_centre_bout_bass"], c["block_bass_upper_bout"], block_top_z, block_bottom_z, step_z, resolution);
	trim_end_block(c["rib_internal_upper_bout"], c["centre_line"], c["block_neck"], 31, 0, 5, 0.1);
	trim_corner_block(c["rib_internal_upper_bout"], c["rib_internal_centre_bout_trbl"], c["block_trbl_upper_bout"], block_top_z, block_bottom_z, step_z, resolution);

}