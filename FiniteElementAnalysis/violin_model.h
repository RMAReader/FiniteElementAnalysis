#ifndef _FEA_VIOLIN_MODEL_H_
#define _FEA_VIOLIN_MODEL_H_


#include <string>
#include <unordered_map>
#include <algorithm>

#include "bspline_utils.h"
#include "geometry.h"

class violin_ribs
{
public:
	std::unordered_map < std::string, geoCURVE2F > curves;
	std::unordered_map < std::string, float > floats;
	std::vector<circle2f> rib_mould_locator_holes;

	void scale_model(double ratio);
	void rotate_model(double angle);
	void translate_model(double angle);
};



class violin_component
{
public:
	violin_component();
	~violin_component();
	violin_component(const violin_component& other);	//copy constructor
	violin_component(violin_component&& other);	//move constructor
	violin_component& operator=(const violin_component& other);	//copy assignment
	violin_component& operator=(violin_component&& other);	//move assignment operator

	std::unordered_map < std::string, geoSURFACE3F > surfaces;
	std::unordered_map < std::string, geoCURVE3F > curves;
	std::unordered_map < std::string, float > floats;

	void rotate_model(double angle);
};

class violin_model
{
public:
	violin_model();
	~violin_model();
	violin_model(const violin_model& other); //copy constructor
	violin_model(violin_model&& other);   //move constructor
	violin_model& operator=(const violin_model& other); //copy assignment
	violin_model& operator=(violin_model&& other);  //move assignment operator


	void scale_model(double ratio);

//private:

	std::string name;
	std::string description;

	violin_ribs ribs;
	violin_component back;
	violin_component belly;

};


#endif _FEA_VIOLIN_MODEL_H_