#include "main.h"


IMPLEMENT_APP(MyApp)

bool MyApp::OnInit()
{

	AppFrame *frame = new AppFrame(wxT("Violin Builder"));

	frame->Show(true);

	return true;
}