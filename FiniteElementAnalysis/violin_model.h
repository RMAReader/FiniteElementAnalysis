#ifndef _FEA_VIOLIN_MODEL_H_
#define _FEA_VIOLIN_MODEL_H_


#include <string>
#include <unordered_map>

#include "bspline_utils.h"


class violin_ribs
{
public:
	std::unordered_map < std::string, CURVE2F* > curves;
	std::unordered_map < std::string, float > floats;
	std::vector<circle2f> rib_mould_locator_holes;

	void scale_model(double ratio);
};

class violin_belly
{
public:
	std::unordered_map < std::string, SURFACE3F* > surfaces;
	std::unordered_map < std::string, CURVE3F* > curves;
	std::unordered_map < std::string, float > floats;
	std::vector<circle2f> rib_mould_locator_holes;

};

class violin_model
{
public:
	violin_model();
	~violin_model();
	

	void scale_model(double ratio);

//private:

	std::string name;
	std::string description;

	violin_ribs* ribs;

};


#endif _FEA_VIOLIN_MODEL_H_