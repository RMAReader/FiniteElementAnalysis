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

	//static void SaveModelXML(const char*, violin_model*, bool);
	static void SaveModelXML(const char*, app_model*, bool);

private:
	
	static std::string IO::ToString(double d);
	static TiXmlElement* IO::NewXmlElement(float x, std::string);
	static TiXmlElement* IO::NewXmlElement(geoVEC2F& point);
	static TiXmlElement* IO::NewXmlElement(geoVEC3F& point);
	static TiXmlElement* IO::NewXmlElement(geoKNOT&, std::string);
	static TiXmlElement* IO::NewXmlElement(geoCURVE2F&, std::string);
	static TiXmlElement* IO::NewXmlElement(geoCURVE3F&, std::string);
	static TiXmlElement* IO::NewXmlElement(geoLATTICE3F&, std::string);
	static TiXmlElement* IO::NewXmlElement(geoSURFACE3F&, std::string);
	static TiXmlElement* IO::NewXmlElement(violin_ribs*, std::string);
	static TiXmlElement* IO::NewXmlElement(violin_component&, std::string);
	static TiXmlElement* IO::NewXmlElement(violin_model&);

	static int IO::NewInteger(TiXmlElement&);
	static float IO::NewFloat(TiXmlElement&);
	static geoVEC2F NewPoint2f(TiXmlElement&);
	static geoVEC3F NewPoint3f(TiXmlElement&);
	static geoKNOT NewKnot(TiXmlElement&);
	static geoCURVE2F NewCurve2f(TiXmlElement&);
	static geoCURVE3F NewCurve3f(TiXmlElement&);
	static geoLATTICE3F NewLattice3f(TiXmlElement&);
	static geoSURFACE3F NewSurface3f(TiXmlElement&);
	static violin_ribs* NewRibs(TiXmlElement&, bool);
	static violin_component NewComponent(TiXmlElement&, bool);
	static violin_model* IO::NewViolin(TiXmlElement*, bool);
	

	static TiXmlElement* IO::NewXmlElement(toolpath_base&);
	static TiXmlElement* IO::NewXmlElement(cnc_tool&);

	static cnc_tool* IO::NewCNCTool(TiXmlElement&, bool);
	static void IO::LoadToolPath(TiXmlElement*, toolpath_base*, bool);
	static toolpath_RibMould* IO::NewToolPath_RibMould(TiXmlElement&, bool);
	static toolpath_TrimBlocksCentreBout* IO::NewToolPath_TrimBlocksCentreBout(TiXmlElement&, bool);
	static toolpath_TrimBlocksEndBouts* IO::NewToolPath_TrimBlocksEndBouts(TiXmlElement&, bool);
	static toolpath_BackRough* IO::NewToolPath_BackRough(TiXmlElement& , bool);
	static toolpath_BackFinish* IO::NewToolPath_BackFinish(TiXmlElement& , bool);
};


#endif _FEA_FILEREADER_H_