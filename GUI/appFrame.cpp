#include "appFrame.h"


AppFrame::AppFrame(const wxString & title)
	: wxFrame(NULL, wxID_ANY, title, wxDefaultPosition, wxSize(600, 600))
{

	wxMenuBar *menubar = new wxMenuBar;
	wxMenu *file = new wxMenu;
	file->Append(wxID_NEW, wxT("&New"));
	file->Append(wxID_OPEN, wxT("&Open"));
	file->Append(wxID_SAVE, wxT("&Save As..."));
	file->Append(wxID_EXIT, wxT("&Exit"));
	menubar->Append(file, wxT("&File"));

	wxMenu *model = new wxMenu;
	model->Append(myID_NEW_CURVE, wxT("&New Curve..."));
	model->Append(myID_DRAG_POINT, wxT("&Drag Point..."));
	menubar->Append(model, wxT("&Model"));

	SetMenuBar(menubar);

	Connect(wxID_NEW, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(AppFrame::OnNew));
	Connect(wxID_OPEN, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(AppFrame::OnOpen));
	Connect(wxID_SAVE, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(AppFrame::OnSave));
	Connect(wxID_EXIT, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(AppFrame::OnExit));
	Connect(myID_NEW_CURVE, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(AppFrame::OnNewCurve));
	Connect(myID_DRAG_POINT, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(AppFrame::OnDragPoint));

	
	CreateStatusBar(3);
	SetStatusText(wxT("Viewpane : (x, y)"), 0);
	SetStatusText(wxT("Pixel : (x, y)"), 1);
	SetStatusText(wxT("Size of line: n"), 2);

	wxToolBar *toolbar = CreateToolBar(wxNO_BORDER | wxHORIZONTAL | wxTB_FLAT, wxID_ANY);
	toolbar->AddRadioTool(myID_NEW_CURVE, wxT("New Curve"), wxBitmap(wxT("C:\\Users\\Lizzie\\Documents\\Richards documents\\pencil18.bmp"), wxBITMAP_TYPE_ANY));
	toolbar->AddRadioTool(myID_DRAG_POINT, wxT("Drag"), wxBitmap(wxT("C:\\Users\\Lizzie\\Documents\\Richards documents\\Drag18.bmp"), wxBITMAP_TYPE_ANY));
	toolbar->Realize();


	canvas = new glCanvas(this, this->GetClientSize());
	const wxGLContext *context = new wxGLContext(canvas);
	bool d = (*canvas).SetCurrent(*context);

	Center();


}

void AppFrame::OnNew(wxCommandEvent& event)
{
	
}
void AppFrame::OnOpen(wxCommandEvent& event)
{

	wxFileDialog * openFileDialog = new wxFileDialog(this, "", "", "", "XML files (*.xml)|*.xml|Dat files (*.dat)|*.dat", wxFD_OPEN, wxDefaultPosition, wxDefaultSize, wxFileDialogNameStr);

	if (openFileDialog->ShowModal() == wxID_OK){
		std::string reverse_path = openFileDialog->GetPath();
		std::reverse(reverse_path.begin(), reverse_path.end());

		std::string extn;
		std::istringstream iss(reverse_path);

		if (std::getline(iss, extn, '.'))
		{
			if (extn == "tad")
			{
				IO::LoadModelJava(openFileDialog->GetPath(), canvas->GetCurves());

				int nc = canvas->GetCurves()->size();
				for (int i = 0; i < nc; i++)
				{
					int np = canvas->GetCurves()->at(i).nPoints();
					for (int j = 0; j < np; j++)
					{
						canvas->GetCurves()->at(i).item(j)->operator*=(0.003f);
					}
				}
			}
			if (extn == "lmx")
			{
				IO::LoadModelXML(openFileDialog->GetPath(), canvas->GetCurves());
			}
		}
		else
		{

		}

	}
}

void AppFrame::OnSave(wxCommandEvent& event)
{

	wxFileDialog * saveFileDialog = new wxFileDialog(this, "","","", "XML files (*.xml)|*.xml|Dat files (*.dat)|*.dat",wxFD_SAVE,wxDefaultPosition,wxDefaultSize,wxFileDialogNameStr);

	if (saveFileDialog->ShowModal() == wxID_OK){
		std::string reverse_path = saveFileDialog->GetPath();
		std::reverse(reverse_path.begin(), reverse_path.end());

		std::string extn;
		std::istringstream iss(reverse_path);
		
		if (std::getline(iss, extn, '.'))
		{
			if (extn == "tad")
			{
				IO::SaveModel(saveFileDialog->GetPath(), *(canvas->GetCurves()));
			}
			if (extn == "lmx")
			{
				IO::SaveModelXML(saveFileDialog->GetPath(), *(canvas->GetCurves()));
			}
		}
		else
		{
			
		}
	}
}

void AppFrame::OnExit(wxCommandEvent& event)
{
	this->Close();
}

void AppFrame::OnNewCurve(wxCommandEvent& event)
{
	canvas->SetModeNewCurve();
}
void AppFrame::OnDragPoint(wxCommandEvent& event)
{
	canvas->SetModeDrag();
}
