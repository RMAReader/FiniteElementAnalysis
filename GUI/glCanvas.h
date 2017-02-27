#ifndef _GLCANVAS_H_
#define _GLCANVAS_H_

#include <wx/wx.h>
#include <GL/glew.h>
#include <wx/glcanvas.h>
#include <vector>
#include <memory>
#include "violin_model.h"



enum mode{
	VIEWING,
	DRAWING_CURVE1,
	DRAWING_CURVE2,
	DRAG_PENDING,
	DRAGGING
};

class glCanvas : public wxGLCanvas
{

private:
	void Render();
	mode currentMode;
	bspline::vec2<int> last_position;
	bspline::vec2<int> leftclick_position;

public:
	glCanvas(wxFrame *parent, const wxSize size);
	~glCanvas();

	void Paintit(wxPaintEvent& event);
	void Sizeit(wxSizeEvent& event);
	void MouseMoved(wxMouseEvent& event);
	void MouseLeftDown(wxMouseEvent& event);
	void MouseLeftUp(wxMouseEvent& event);
	void MouseRightClick(wxMouseEvent& event);
	void SetModeDrag();
	void SetModeNewCurve();
	void AddCurve(const curve2f& curve);
	std::vector<curve2f>* GetCurves();
	//vector of curves
	std::vector<curve2f>* curve;
	std::vector<offsetCurve2f>* offsetCurve;

	//vector of pointers to selected vertices
	std::vector<point2f*>* selection;

	float clip_minx;
	float clip_miny;
	float clip_maxx;
	float clip_maxy;

	void SetClipping();

};



#endif _GLCANVAS_H_