#include "FileReader.h"

void IO::ReadInput(const char* filepath, BSplineSolid* solid, double* D, double& density)
{
	
		std::ifstream input;
		input.open(filepath, std::ios::in);
		if (input.is_open()){
			std::string line;
			std::string::size_type sz;
			size_t start, l;
			while (std::getline(input, line)){

				if (line.find("Input") != 0) continue;
				std::getline(input, line);
				if (line.find("Density(kgm3)") == 0){
					start = line.find(" ");
					l = line.length() - start;
					density = stod(line.substr(start, l));
				}
				std::getline(input, line);
				std::getline(input, line);

				int counter = 0;
				while (counter < 36){
					do{
						std::getline(input, line, ' ');
					} while (line.length() == 0 || line =="\n");
					D[counter] = stod(line);
					counter++;
				}
				std::getline(input, line);

				int p;
				std::getline(input, line, ' ');
				if (line.find("p:") == 0){
					std::getline(input, line);
					p = stoi(line);
				}
				int q;
				std::getline(input, line, ' ');
				if (line.find("q:") == 0){
					std::getline(input, line);
					q = stoi(line);
				}
				int r;
				std::getline(input, line, ' ');
				if (line.find("r:") == 0){
					std::getline(input, line);
					r = stoi(line);
				}
				int lknotx;
				std::getline(input, line, ' ');
				if (line.find("lknotx:") == 0){
					std::getline(input, line);
					lknotx = stoi(line);
				}
				int lknoty;
				std::getline(input, line, ' ');
				if (line.find("lknoty:") == 0){
					std::getline(input, line);
					lknoty = stoi(line);
				}
				int lknotz;
				std::getline(input, line, ' ');
				if (line.find("lknotz:") == 0){
					std::getline(input, line);
					lknotz = stoi(line);
				}
				int nx;
				std::getline(input, line, ' ');
				if (line.find("nx:") == 0){
					std::getline(input, line);
					nx = stoi(line);
				}
				int ny;
				std::getline(input, line, ' ');
				if (line.find("ny:") == 0){
					std::getline(input, line);
					ny = stoi(line);
				}
				int nz;
				std::getline(input, line, ' ');
				if (line.find("nz:") == 0){
					std::getline(input, line);
					nz = stoi(line);
				}


				solid->initialise(p, q, r, nx, ny, nz);
				for (int i = 0; i < lknotx; i++){
					std::getline(input, line);
					start = line.find(" ");
					l = line.length() - start;
					solid->knotx[i] = stod(line.substr(start, l));
				}
				for (int i = 0; i < lknoty; i++){
					std::getline(input, line);
					start = line.find(" ");
					l = line.length() - start;
					solid->knoty[i] = stod(line.substr(start, l));
				}
				for (int i = 0; i < lknotz; i++){
					std::getline(input, line);
					start = line.find(" ");
					l = line.length() - start;
					solid->knotz[i] = stod(line.substr(start, l));
				}

				for (int i = 0; i < nx*ny*nz; i++){
					std::getline(input, line);
					solid->cx[i] = stod(line, &sz);
					solid->cy[i] = stod(line.substr(sz, sz));
					solid->cz[i] = stod(line.substr(2 * sz, sz));
				}

			}
			input.close();
		}
		else{
			std::cout << "Unable to read input file";
		}
	
}



bool IO::LoadModelJava(const char* filepath, std::vector<geometry::bspline::curve<geometry::vector<float, 2>, double>>& curve, std::vector<geometry::bspline::surface<geometry::vector<float, 3>, double,double>>& surfaces, bool verbose)
{
	std::ifstream input;
	input.open(filepath, std::ios::in);
	if (input.is_open()){
		std::string line, word;


		while (std::getline(input, line, '\n')){

			if (line.find("BeginGraphicBSplineSurface") == 0)
			{
				int m, n;
				geometry::bspline::knot<double> knotx;
				geometry::bspline::knot<double> knoty;
				std::vector<float> Px;
				std::vector<float> Py;
				std::vector<float> Pz;
				geoLATTICE3F lattice;
				m = 0;
				n = 0;
				while (line.find("End") != 0)
				{
					std::getline(input, line, '\n');
					std::istringstream iss(line);

					if (line.find("p") == 0){
						std::getline(iss, word, ':');
						std::getline(iss, word, '\n');
						knotx.p = stoi(word);
					}
					else if (line.find("q") == 0){
						std::getline(iss, word, ':');
						std::getline(iss, word, '\n');
						knoty.p = stoi(word);
					}
					else if (line.find("knotx") == 0)
					{
						std::getline(iss, word, ':');
						while (std::getline(iss, word, ','))
						{
							knotx._data.push_back(stod(word));
							if (verbose){ std::cout << "knotx[" << knotx._data.size() << "] = " << word << std::endl; }
						}
					}
					else if (line.find("knoty") == 0)
					{
						std::getline(iss, word, ':');
						while (std::getline(iss, word, ','))
						{
							knoty._data.push_back(stod(word));
							if (verbose){ std::cout << "knoty[" << knoty._data.size() << "] = " << word << std::endl; }
						}
					}
					else if (line.find("Px[") == 0)
					{
						m++;
						std::getline(iss, word, ':');
						while (std::getline(iss, word, ','))
						{
							Px.push_back((float)stod(word));
							if (verbose){ std::cout << "Px[" << Px.size() << "] = " << word << std::endl; }
						}
					}
					else if (line.find("Py[") == 0)
					{
						std::getline(iss, word, ':');
						while (std::getline(iss, word, ','))
						{
							Py.push_back((float)stod(word));
							if (verbose){ std::cout << "Py[" << Py.size() << "] = " << word << std::endl; }
						}
					}
					else if (line.find("Pz[") == 0)
					{
						std::getline(iss, word, ':');
						while (std::getline(iss, word, ','))
						{
							Pz.push_back((float)stod(word));
							if (verbose){ std::cout << "Pz[" << Pz.size() << "] = " << word << std::endl; }
						}
					}
				}
				if (Px.size() == Py.size() && Px.size() == Pz.size())
				{
					n = Px.size() / m;
					geoLATTICE3F lattice(m, n);
					for (int i = 0; i < Px.size(); i++){
						lattice.data[i] = geoVEC3F(std::array < float, 3 > {{Px[i], Py[i], Pz[i]}});
					}
					surfaces.push_back(geoSURFACE3F(knotx.p,knoty.p, knotx, knoty, lattice));
				}

			}
			else if (line.find("BeginGraphicBSplineCurve:") == 0)
			{
				geometry::bspline::curve<geometry::vector<float,2>, double>  new_curve;
				
				while (line.find("End") != 0){
					
					std::getline(input, line, '\n');
					std::istringstream iss(line);

					if (line.find("p") == 0){
						std::getline(iss, word, ':');
						std::getline(iss, word, '\n');
						new_curve._p = stoi(word);
					}
					else if (line.find("knot") == 0)
					{
						std::getline(iss, word, ':');
						while (std::getline(iss, word, ','))
						{
							new_curve._knot._data.push_back(stod(word));
							if (verbose){ std::cout << "knot[" << new_curve._knot._data.size() << "] = " << word << std::endl; }
						}
					}
					else if (line.find("P[") == 0)
					{
						geometry::vector<float,2> p;
						std::getline(iss, word, ':');
						std::getline(iss, word, ',');
						p[0] = stod(word);
						std::getline(iss, word, ',');
						p[1] = stod(word);
						/*std::getline(iss, word, ',');
						z = stod(word);*/
						new_curve._points.push_back(p);
					}
				}

				//point2f* points_final = new point2f[(int)points.size()];
				//double* knot_final = new 

				curve.push_back(new_curve);
			}
			 
		}
		input.close();
		return true;
	}
	else{
		return false;
	}
}




void IO::SaveModelXML(const char* filepath, violin_model* violin, bool verbose){
	
	TiXmlDocument doc;
	TiXmlDeclaration * decl = new TiXmlDeclaration("1.0", "", "");
	doc.LinkEndChild(decl);

	TiXmlElement* violin_element = new TiXmlElement("Violin");
	doc.LinkEndChild(violin_element);

	violin_element->SetAttribute("name", violin->name);
	TiXmlElement* description = new TiXmlElement("Description");
	TiXmlText* text = new TiXmlText(violin->description);
	description->LinkEndChild(text);
	violin_element->LinkEndChild(description);

	std::string name = "violin_ribs";
	TiXmlElement* ribs = IO::NewXmlElement(&(violin->ribs), name);
	violin_element->LinkEndChild(ribs);

	doc.SaveFile(filepath);

}


std::string IO::ToString(double d)
{
	std::ostringstream s; s << d; return s.str();
}

TiXmlElement* IO::NewXmlElement(float x, std::string name){
	TiXmlElement* result = new TiXmlElement("float");
	result->SetAttribute("name", name);
	result->SetDoubleAttribute("value", x);
	return result;
}



TiXmlElement* IO::NewXmlElement(geoVEC2F& point)
{
	TiXmlElement* result = new TiXmlElement("vector2f");
	result->SetDoubleAttribute("0", point[0]);
	result->SetDoubleAttribute("1", point[1]);
	return result;
}



TiXmlElement* IO::NewXmlElement(geoCURVE2F& curve, std::string name)
{
	TiXmlElement* result = new TiXmlElement("curve2f");
	result->SetAttribute("name", name);

	TiXmlElement* order = new TiXmlElement("order");
	order->SetAttribute("value", ToString(curve._p));
	result->LinkEndChild(order);

	TiXmlElement* knot = new TiXmlElement("knot");
	result->LinkEndChild(knot);
	for (int i = 0; i < curve._knot.size(); i++)
	{
		TiXmlElement* value = new TiXmlElement("double");
		value->SetDoubleAttribute("value", curve._knot[i]);
		knot->LinkEndChild(value);
	}

	TiXmlElement* points = new TiXmlElement("points");
	result->LinkEndChild(points);
	for (int i = 0; i < curve._points.size(); i++)
	{
		TiXmlElement* value = NewXmlElement(curve._points[i]);
		points->LinkEndChild(value);
	}
	return result;
}


TiXmlElement* IO::NewXmlElement(violin_ribs* ribs, std::string name)
{
	TiXmlElement* result = new TiXmlElement("violin_ribs");
	for (auto it = ribs->curves.begin(); it != ribs->curves.end(); ++it){
		TiXmlElement* e = NewXmlElement(it->second, it->first);
		result->LinkEndChild(e);
	}
	for (auto it = ribs->floats.begin(); it != ribs->floats.end(); ++it){
		result->LinkEndChild(NewXmlElement(it->second, it->first));
	}
	return result;
}




app_model* IO::LoadModelXML(const char* filepath,  bool verbose)
{
	TiXmlDocument doc(filepath);
	if (!doc.LoadFile()) return false;

	app_model* model = new app_model();

	for (TiXmlElement* e = doc.FirstChildElement(); e != NULL; e = e->NextSiblingElement())
	{
		std::string type = e->Value();
		if (type == "violin_model")
		{
			model->violin = NewViolin(e, verbose);
		}
		else if (type=="toolpath")
		{
			std::string subtype = e->Attribute("type");
			if (subtype == "toolpath_RibMould")
			{
				toolpath_RibMould* path = NewToolPath_RibMould(*e, verbose);
				model->paths.insert(std::make_pair(subtype, path));
			}
			else if (subtype == "toolpath_TrimBlocksEndBouts")
			{
				toolpath_TrimBlocksEndBouts* path = NewToolPath_TrimBlocksEndBouts(*e, verbose);
				model->paths.insert(std::make_pair(subtype, path));
			}
			else if (subtype == "toolpath_TrimBlocksCentreBout")
			{
				toolpath_TrimBlocksCentreBout* path = NewToolPath_TrimBlocksCentreBout(*e, verbose);
				model->paths.insert(std::make_pair(subtype, path));
			}
			else if (subtype == "toolpath_BackRough")
			{
				toolpath_BackRough* path = NewToolPath_BackRough(*e, verbose);
				model->paths.insert(std::make_pair(subtype, path));
			}
		}
	}

	//link each toolpath to the violin object
	for (std::pair<std::string,toolpath_base*> path: model->paths){
		path.second->violin = model->violin;
	}
	return model;
}




violin_model* IO::NewViolin(TiXmlElement* doc, bool verbose)
{
	violin_model* violin = new violin_model();
	violin->name = doc->Attribute("name");

	for (TiXmlElement* e = doc->FirstChildElement(); e != NULL; e = e->NextSiblingElement())
	{
		std::string type = e->Value();
		if (type == "violin_ribs")
		{
			violin->ribs = *NewRibs(*e, verbose);
		}
		else if (type == "description"){
			violin->description = e->GetText();
		}
	}
	return violin;
}


violin_ribs* IO::NewRibs(TiXmlElement& root,  bool verbose)
{
	violin_ribs* ribs = new violin_ribs();
	
	for (TiXmlElement* e = root.FirstChildElement(); e != NULL; e = e->NextSiblingElement())
	{
		std::string type = e->Value();
		if (type == "curve2f")
		{
			std::string name = e->Attribute("name");
			ribs->curves.insert(std::make_pair(name, NewCurve2f(*e)));

		}
		else if (type == "float"){

			std::string name = e->Attribute("name");
			double d;
			bool invalid = e->QueryValueAttribute("value", &d);

			ribs->floats.insert(std::make_pair(name, d));
		}
	}
	return ribs;
}



int IO::NewInteger(TiXmlElement& e)
{
	std::string name = e.Value();

	int d;
	bool invalid = e.QueryIntAttribute("value", &d);
	if (invalid){ return 0; }
	else { return d; }
}
float IO::NewFloat(TiXmlElement& e)
{
	std::string name = e.Value();

	double d;
	bool invalid = e.QueryValueAttribute("value", &d);
	if (invalid){ return 0; }
	else { return d; }
}

//VEC2F IO::NewPoint2f(TiXmlElement& e)
//{
//	point2f point;
//	e.QueryValueAttribute("x", &(point.x));
//	e.QueryValueAttribute("y", &(point.y));
//	return point;
//}
//CURVE2F* IO::NewCurve2f(TiXmlElement& root)
//{
//	int p;
//	double k;
//	std::vector<double> knot;
//	std::vector<bspline::vec2<float>> points;
//	std::string name;
//
//	for (TiXmlElement* e = root.FirstChildElement(); e != NULL; e = e->NextSiblingElement())
//	{
//		name = e->Value();
//		if (name == "order")
//		{
//			bool valid = e->QueryValueAttribute("value", &p);
//		}
//		else if (name == "knot")
//		{
//			for (TiXmlElement* ke = e->FirstChildElement(); ke != NULL; ke = ke->NextSiblingElement())
//			{
//				
//				double d;
//				bool valid = ke->QueryValueAttribute("value", &d);
//				knot.push_back(d);
//			}
//		}
//		else if (name == "points")
//		{
//			for (TiXmlElement* ke = e->FirstChildElement(); ke != NULL; ke = ke->NextSiblingElement())
//			{
//				points.push_back(NewPoint2f(*ke));
//			}
//		}
//	}
//	return new CURVE2F(p, (int)points.size(), (int)knot.size(), knot.data(), points.data());
//}

geoVEC2F IO::NewPoint2f(TiXmlElement& e)
{
	geoVEC2F point;
	e.QueryValueAttribute("x", &(point[0]));
	e.QueryValueAttribute("y", &(point[1]));
	return point;
}
geoCURVE2F IO::NewCurve2f(TiXmlElement& root)
{
	geometry::bspline::knot<double> knot;
	std::vector<geometry::vector<float,2>> points;
	std::string name;

	for (TiXmlElement* e = root.FirstChildElement(); e != NULL; e = e->NextSiblingElement())
	{
		name = e->Value();
		if (name == "order")
		{
			bool valid = e->QueryValueAttribute("value", &(knot.p));
		}
		else if (name == "knot")
		{
			for (TiXmlElement* ke = e->FirstChildElement(); ke != NULL; ke = ke->NextSiblingElement())
			{

				double d;
				bool valid = ke->QueryValueAttribute("value", &d);
				knot._data.push_back(d);
			}
			knot.normalise();
		}
		else if (name == "points")
		{
			for (TiXmlElement* ke = e->FirstChildElement(); ke != NULL; ke = ke->NextSiblingElement())
			{
				points.push_back(NewPoint2f(*ke));
			}
		}
	}
	return geoCURVE2F(knot.p, knot, points);
}



cnc_tool* IO::NewCNCTool(TiXmlElement& e, bool verbose){

	std::string type;
	std::string name;
	
	cnc_tool* tool = new cnc_tool();

	tool->name = e.Attribute("name");
	for (TiXmlElement* f = e.FirstChildElement(); f != NULL; f = f->NextSiblingElement())
	{
		type = f->Value();
		if (type == "float")
		{
			name = f->Attribute("name");
			if (name == "diameter"){
				tool->diameter = NewFloat(*f);
			}
		}
		else if (type == "string"){
			tool->type = f->Attribute("type");
		}
	}
	return tool;
}


void IO::LoadToolPath(TiXmlElement* root,toolpath_base* toolpath, bool verbose)
{
	for (TiXmlElement* e = root->FirstChildElement(); e != NULL; e = e->NextSiblingElement())
	{
		std::string type;
		std::string name;
		type = e->Value();
		if (type == "string")
		{
			name = e->Attribute("name");
			if (name == "gcodefile"){ toolpath->gcode_filepath = e->Attribute("value"); }
		}
		else if (type == "tool")
		{
			toolpath->tool = NewCNCTool(*e, verbose);
		}
		else if (type == "float")
		{
			name = e->Attribute("name");
			float value = NewFloat(*e);
			toolpath->parameters.insert(std::make_pair(name, value));
		}
	}
}


toolpath_RibMould* IO::NewToolPath_RibMould(TiXmlElement& root, bool verbose)
{
	toolpath_RibMould* toolpath = new toolpath_RibMould();
	IO::LoadToolPath(&root, toolpath, verbose);
	return toolpath;
}
toolpath_TrimBlocksCentreBout* IO::NewToolPath_TrimBlocksCentreBout(TiXmlElement& root, bool verbose)
{
	toolpath_TrimBlocksCentreBout* toolpath = new toolpath_TrimBlocksCentreBout();
	IO::LoadToolPath(&root, toolpath, verbose);
	return toolpath;
}
toolpath_TrimBlocksEndBouts* IO::NewToolPath_TrimBlocksEndBouts(TiXmlElement& root, bool verbose)
{
	toolpath_TrimBlocksEndBouts* toolpath = new toolpath_TrimBlocksEndBouts();
	IO::LoadToolPath(&root, toolpath, verbose);
	return toolpath;
}
toolpath_BackRough* IO::NewToolPath_BackRough(TiXmlElement& root, bool verbose)
{
	toolpath_BackRough* toolpath = new toolpath_BackRough();
	IO::LoadToolPath(&root, toolpath, verbose);
	return toolpath;
}
