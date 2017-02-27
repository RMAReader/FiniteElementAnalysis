#pragma once
#include <string>


#ifdef FEADLL_EXPORTS
#define FEADLL_API __declspec(dllexport) 
#else
#define FEADLL_API __declspec(dllimport) 
#endif

class HarwellBoeing
{
public:
	HarwellBoeing();
	~HarwellBoeing();

	static int Save(std::string filename, int NROW, int NCOL, int NNZERO, double VALUES[], int ROWIND[], int COLPTR[]);
	static int Load(std::string filename, double Matrix[]);
};

