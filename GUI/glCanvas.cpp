#include "glCanvas.h"




glCanvas::glCanvas(wxFrame *parent, const wxSize size)
	:wxGLCanvas(parent, wxID_ANY, 0, wxDefaultPosition, size, 0, wxT("GLCanvas"), wxNullPalette){

	Bind(wxEVT_PAINT, &glCanvas::Paintit, this);
	Bind(wxEVT_SIZE, &glCanvas::Sizeit, this);
	Bind(wxEVT_MOTION, &glCanvas::MouseMoved, this);
	Bind(wxEVT_LEFT_UP, &glCanvas::MouseLeftUp, this);
	Bind(wxEVT_LEFT_DOWN, &glCanvas::MouseLeftDown, this);
	Bind(wxEVT_RIGHT_UP, &glCanvas::MouseRightClick, this);

	this->SetCursor(wxCursor(wxCURSOR_CROSS));
	
	currentMode = VIEWING;
	selection = new std::vector<point2f*>();
	curve = new std::vector<curve2f>();
	offsetCurve = new std::vector<offsetCurve2f>();

	// Initialize GLEW
	//if (glewInit() != GLEW_OK) {
	//	fprintf(stderr, "Failed to initialize GLEW\n");
	//	this->GetParent()->Close();
	//}

	std::ofstream log;
	log.open("glCanvas Log.txt");
	log.close();


}


glCanvas::~glCanvas()
{
}


void glCanvas::Paintit(wxPaintEvent& WXUNUSED(event)){
	Render();
}
void glCanvas::Sizeit(wxSizeEvent& WXUNUSED(event)){
	this->Refresh();
}
void glCanvas::MouseMoved(wxMouseEvent& event){
	
	wxFrame *myParent = (wxFrame *) this->GetParent();
	int pixel_x = GetSize().GetWidth();
	int pixel_y = GetSize().GetHeight();

	float x = clip_minx * (1.0 - (float)event.GetX() / pixel_x) + clip_maxx * (float) event.GetX() / pixel_x;
	float y = clip_miny * (float) event.GetY() / pixel_y + clip_maxy * (1.0 - (float)event.GetY() / pixel_y);

	float dx = (event.GetX() - last_position.x) * (clip_maxx - clip_minx) / GetSize().GetWidth();
	float dy = -(event.GetY() - last_position.y) * (clip_maxy - clip_miny) / GetSize().GetHeight();

	if (currentMode == DRAWING_CURVE1 || currentMode == DRAWING_CURVE2 || currentMode == DRAGGING)
	{
		for (int i = 0; i < selection->size(); i++)
		{ 
			((*selection)[i])->Translate(dx,dy);
		}
	}

	last_position = bspline::vec2<int>(event.GetX(), event.GetY());
	

	this->Refresh();

	myParent->SetStatusText(wxString::Format("Viewpane: (%0f,%0f)", x, y), 0);
	myParent->SetStatusText(wxString::Format("Pixel: (%d,%d)", pixel_x, pixel_y), 1);

}

void glCanvas::MouseLeftDown(wxMouseEvent& event)
{
	int pixel_x = GetSize().GetWidth();
	int pixel_y = GetSize().GetHeight();

	float x = clip_minx * (1.0 - (float)event.GetX() / pixel_x) + clip_maxx * (float)event.GetX() / pixel_x;
	float y = clip_miny * (float)event.GetY() / pixel_y + clip_maxy * (1.0 - (float)event.GetY() / pixel_y);
	
	if (currentMode == DRAG_PENDING)
	{
		point2f cursor(x, y);
		float d = 0.05f;
		for (int i = 0; i < curve->size(); i++)
		{
			for (int j = 0; j < (*curve)[i].nPoints(); j++)
			{
				float e = cursor.GetDistance((*curve)[i].get(j));
				if (e < d)
				{
					d = e;
					selection->clear();
					selection->push_back((*curve)[i].item(j));
					currentMode = DRAGGING;
				}
			}
		}
	}

	this->Refresh();

}

void glCanvas::MouseLeftUp(wxMouseEvent& event)
{
	wxFrame *myParent = (wxFrame *) this->GetParent();
	int pixel_x = GetSize().GetWidth();
	int pixel_y = GetSize().GetHeight();

	float x = clip_minx * (1.0 - (float)event.GetX() / pixel_x) + clip_maxx * (float)event.GetX() / pixel_x;
	float y = clip_miny * (float)event.GetY() / pixel_y + clip_maxy * (1.0 - (float)event.GetY() / pixel_y);

	last_position = bspline::vec2<int>(event.GetX(), event.GetY());

	if (currentMode == VIEWING)
	{
		curve->push_back(curve2f(2, point2f(x, y)));
		curve->back().append(point2f(x, y));
		curve->back().append(point2f(x, y));

		selection->push_back(curve->back().item(2));
		selection->push_back(curve->back().item(1));
		currentMode = DRAWING_CURVE1;

		//temp code to add offset curve as a test
		//offsetCurve->push_back(offsetCurve2f(&(curve->back()), 0.05));
	}
	else if (currentMode == DRAWING_CURVE1)
	{
		selection->pop_back();
		currentMode = DRAWING_CURVE2;
	}
	else if (currentMode == DRAWING_CURVE2)
	{
		curve2f *currentCurve = (curve2f*)selection->front()->GetParent();
		currentCurve->append(point2f(x, y));
		int n = currentCurve->nPoints();
		selection->pop_back();
		selection->push_back(currentCurve->item(n - 1));
	}
	else if (currentMode == DRAGGING)
	{
		selection->clear();
		currentMode = DRAG_PENDING;
	}

	this->Refresh();

}


void glCanvas::MouseRightClick(wxMouseEvent& event)
{
	if (currentMode == DRAWING_CURVE2)
	{
		selection->clear();
		currentMode = VIEWING;
	}
	
}


void glCanvas::Render()
{
	//wxPaintDC(this);
	glClearColor(1.0, 1.0, 1.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT);
	
	//set viewport and clipping area
	SetClipping();

	//glClear(GL_COLOR_BUFFER_BIT);   // Clear the color buffer with current clearing color


	for (int c = 0; c < curve->size(); c++)
	{
		glBegin(GL_LINE_STRIP);            
		glColor3f(0.3f, 0.3f, 0.3f); // Gray
		for (int i = 0; i < (*curve)[c].nPoints(); i++){
			glVertex2f((GLfloat)(*curve)[c].get(i).x, (GLfloat)(*curve)[c].get(i).y);
		}
		glEnd();

		glBegin(GL_LINE_STRIP);            
		glColor3f(0.0f, 0.0f, 1.0f); // Blue
		for (int i = 0; i <= 1000; i++){
			double t = ((*curve)[c].minParam() * (1000 - i) + (*curve)[c].maxParam() * i) / 1000;
			point2f p = (*curve)[c].evaluate(t);
			glVertex2f((GLfloat)p.x, (GLfloat)p.y);
		}
		glEnd();
	}

	for (int c = 0; c < offsetCurve->size(); c++)
	{
		glBegin(GL_LINE_STRIP);
		glColor3f(1.0f, 0.0f, 0.0f); // Blue
		offsetCurve2f* current = &(*offsetCurve)[c];
		for (int i = 0; i <= 1000; i++){
			double t = (current->GetBase()->minParam() * (1000 - i) + current->GetBase()->maxParam() * i) / 1000;
			point2f p = current->evaluate(t);
			glVertex2f((GLfloat)p.x, (GLfloat)p.y);
		}
		glEnd();
	}


	glFlush();  // Render now
	
	SwapBuffers();
}

void glCanvas::SetClipping()
{
	glViewport(0, 0, (GLint)GetSize().GetWidth(), (GLint)GetSize().GetHeight());
	glMatrixMode(GL_PROJECTION);      // Select the Projection matrix for operation
	glLoadIdentity();
	if (GetSize().GetHeight() < GetSize().GetWidth())
	{
		float aspectRatio = (float)GetSize().GetHeight() / (float)GetSize().GetWidth();
		clip_minx = -1.0f;
		clip_maxx = 1.0f;
		clip_miny = -aspectRatio;
		clip_maxy = aspectRatio;

		gluOrtho2D(clip_minx, clip_maxx, clip_miny, clip_maxy);
	}
	else
	{
		float aspectRatio = (float)GetSize().GetWidth() / (float)GetSize().GetHeight();
		clip_minx = -aspectRatio;
		clip_maxx = aspectRatio;
		clip_miny = -1.0;
		clip_maxy = 1.0;

		gluOrtho2D(clip_minx, clip_maxx, clip_miny, clip_maxy);
	}
}


void glCanvas::SetModeDrag()
{
	selection->clear();
	currentMode = DRAG_PENDING;
}
void glCanvas::SetModeNewCurve()
{
	selection->clear();
	currentMode = VIEWING;
}
void glCanvas::AddCurve(const curve2f& curve)
{
	this->curve->push_back(curve);
}
std::vector<curve2f>* glCanvas::GetCurves()
{
	return curve;
}