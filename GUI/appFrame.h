#ifndef _APPFRAME_H_
#define _APPFRAME_H_

#include <wx/wx.h>
#include <wx/frame.h>
#include "glCanvas.h"
#include <iostream>
#include <fstream>
#include <algorithm>

#include "FileReader.h"

enum MyID{
	myID_NEW_CURVE,
	myID_DRAG_POINT
};


class AppFrame : public wxFrame
{
private:
	glCanvas* canvas;

public:
	AppFrame(const wxString& title);

	void OnNew(wxCommandEvent& event);
	void OnOpen(wxCommandEvent& event);
	void OnSave(wxCommandEvent& event);
	void OnExit(wxCommandEvent& event);
	void OnNewCurve(wxCommandEvent& event);
	void OnDragPoint(wxCommandEvent& event);



};

#endif _APPFRAME_H_