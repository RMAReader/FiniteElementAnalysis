#ifndef _APP_MODEL_H_
#define _APP_MODEL_H_

#include <unordered_map>

#include "violin_model.h"
#include "toolpath_base.h"

class app_model
{
public:
	
	violin_model* violin = nullptr;
	std::unordered_map < std::string, toolpath_base* > paths;


};

#endif _APP_MODEL_H_
