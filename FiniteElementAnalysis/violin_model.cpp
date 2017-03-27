#include "violin_model.h"


violin_model::violin_model()
{
}


violin_model::~violin_model()
{
}

//copy constructor
violin_model::violin_model(const violin_model& other)
{
	*this = other;
}

//copy assignment
violin_model& violin_model::operator=(const violin_model& other)
{
	if (this != &other)
	{
		//call delete[] on any pointer class members to free existing resource
		name = other.name;
		description = other.description;
		ribs = other.ribs;
		back = other.back;
		belly = other.belly;

	}
	return *this;
}

//move constructor
violin_model::violin_model(violin_model&& other)
{
	*this = std::move(other);
}

//move assignment operator
violin_model& violin_model::operator=(violin_model&& other)
{
	if (this != &other)
	{
		name = std::move(other.name);
		description = std::move(other.description);
		ribs = std::move(other.ribs);
		back = std::move(other.back);
		belly = std::move(other.belly);
	}
	return *this;
}






void violin_ribs::scale_model(double ratio){

	//for (auto it = curves.begin(); it != curves.end(); ++it){
	//	for (int i = 0; i < it->second._points.size(); i++){
	//		it->second._points[i]*=ratio;
	//	}
	//}
	for (auto it = floats.begin(); it != floats.end(); ++it){
		it->second *=ratio;
	}
}
void violin_model::scale_model(double ratio){
	ribs.scale_model(ratio);
}





