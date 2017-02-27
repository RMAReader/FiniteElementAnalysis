#pragma once

#ifdef FEADLL_EXPORTS
#define FEADLL_API __declspec(dllexport) 
#else
#define FEADLL_API __declspec(dllimport) 
#endif


class Matrix
{
public:
	Matrix();
	~Matrix();
};

