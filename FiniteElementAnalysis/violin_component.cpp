#include "violin_model.h"


violin_component::violin_component()
{
}


violin_component::~violin_component()
{
}

//copy constructor
violin_component::violin_component(const violin_component& other)
{
	*this = other;
}

//copy assignment
violin_component& violin_component::operator=(const violin_component& other)
{
	if (this != &other)
	{
		//call delete[] on any pointer class members to free existing resource
		curves = other.curves;
		surfaces = other.surfaces;
		floats = other.floats;
	}
	return *this;
}

//move constructor
violin_component::violin_component(violin_component&& other)
{
	*this = std::move(other);
}

//move assignment operator
violin_component& violin_component::operator=(violin_component&& other)
{
	if (this != &other)
	{
		curves = std::move(other.curves);
		surfaces = std::move(other.surfaces);
		floats = std::move(other.floats);
	}
	return *this;
}


void violin_component::rotate_model(double angle){

	double m1 = sin(angle);
	double m2 = cos(angle);
	for (auto item : surfaces)
	{
		for (int i = 0; i < item.second->points.data.size(); i++){
			VEC3F& p = item.second->points.data[i];
			p = VEC3F(m2*p.x + m1*p.y, m1*p.x - m2*p.y, p.z);
		}
	}
}