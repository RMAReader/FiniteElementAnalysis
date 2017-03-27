#ifndef _FEA_FILEREADER_H_
#define _FEA_FILEREADER_H_

#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <exception>

#include "tinyxml.h"
#include "BSplineSolid.h"
#include "violin_model.h"
#include "toolpath_base.h"
#include "app_model.h"
#include "geometry.h"


class IO
{
public:
	
	static void ReadInput(const char*, BSplineSolid*, double*, double&);
	static bool LoadModelJava(const char*, std::vector<geoCURVE2F>&, std::vector<geoSURFACE3F>&, bool);
	static app_model* LoadModelXML(const char*, bool);

	//static void SaveModel(const char*, std::vector<geoCURVE2F>&);
	//static void SaveModelXML(const char*, std::vector<geoCURVE2F>&);
	static void SaveModelXML(const char*, violin_model*, bool);

private:
	
	static std::string IO::ToString(double d);
	static TiXmlElement* IO::NewXmlElement(float x, std::string);
	static TiXmlElement* IO::NewXmlElement(geoVEC2F& point);
	static TiXmlElement* IO::NewXmlElement(geoCURVE2F&, std::string);
	static TiXmlElement* IO::NewXmlElement(violin_ribs*,std::string);
	
	static int IO::NewInteger(TiXmlElement&);
	static float IO::NewFloat(TiXmlElement&);
	static geoVEC2F NewPoint2f(TiXmlElement&);
	static geoCURVE2F NewCurve2f(TiXmlElement&);
	static violin_ribs* NewRibs(TiXmlElement&, bool);
	static violin_model* IO::NewViolin(TiXmlElement*, bool);
	static cnc_tool* IO::NewCNCTool(TiXmlElement&, bool);

	static void IO::LoadToolPath(TiXmlElement*, toolpath_base*, bool);
	static toolpath_RibMould* IO::NewToolPath_RibMould(TiXmlElement&, bool);
	static toolpath_TrimBlocksCentreBout* IO::NewToolPath_TrimBlocksCentreBout(TiXmlElement&, bool);
	static toolpath_TrimBlocksEndBouts* IO::NewToolPath_TrimBlocksEndBouts(TiXmlElement&, bool);
	static toolpath_BackRough* IO::NewToolPath_BackRough(TiXmlElement& , bool);
};


#endif _FEA_FILEREADER_H_