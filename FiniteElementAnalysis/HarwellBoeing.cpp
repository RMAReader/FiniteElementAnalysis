#include "stdafx.h"
#include "HarwellBoeing.h"
#include <string>
#include "stdio.h"
#include "iostream"
#include <fstream>


using namespace std;

HarwellBoeing::HarwellBoeing(){}

HarwellBoeing::~HarwellBoeing(){}

int HarwellBoeing::Save(string filename, int NROW, int NCOL, int NNZERO, double VALUES[], int ROWIND[], int COLPTR[]){
	string line;
	char* TITLE = "Matrix stored in Harwell-Boeing format";
	string KEY = "Sym";
	int TOTCRD;
	int PTRCRD;
	int INDCRD;
	int VALCRD;
	int RHSCRD = 0;
	string MXTYPE = "RSA";
	int NELTVL = 0;
	string PTRFMT = "(1I9)";
	string INDFMT = "(1I9)";
	string VALFMT = "(1F15.8)";
	string RHSFMT = "(1F15.8)";

	PTRCRD = NCOL;
	INDCRD = NNZERO;
	VALCRD = NNZERO;
	TOTCRD = PTRCRD + INDCRD + VALCRD;

	


	FILE * pFile;
	
	char* fileNameChar = const_cast<char*>(filename.c_str());
	try {
		fopen_s(&pFile,fileNameChar, "w");

		//write header block
		//line 1
		fprintf(pFile, "%72c", TITLE);
		fprintf(pFile, "%8c", KEY);
		fprintf(pFile, "\r\n");
		//line 2
		fprintf(pFile, "%14d", TOTCRD);
		fprintf(pFile, "%14d", PTRCRD);
		fprintf(pFile, "%14d", INDCRD);
		fprintf(pFile, "%14d", VALCRD);
		fprintf(pFile, "%14d", RHSCRD);
		fprintf(pFile, "\r\n");
		//line 3
		fprintf(pFile, "%3c", MXTYPE);
		fprintf(pFile, "           ");
		fprintf(pFile, "%14d", NROW);
		fprintf(pFile, "%14d", NCOL);
		fprintf(pFile, "%14d", NNZERO);
		fprintf(pFile, "%14d", NELTVL);
		fprintf(pFile, "\r\n");
		//line 4
		fprintf(pFile, "%16s", PTRFMT);
		fprintf(pFile, "%16s", INDFMT);
		fprintf(pFile, "%20s", VALFMT);
		fprintf(pFile, "%20s", RHSFMT);
		fprintf(pFile, "\r\n");

		//block of col pointers
		for (int i = 0; i<NCOL; i++) fprintf(pFile, "%010d\r\n", COLPTR[i]);
		//block of row indices
		for (int i = 0; i<NROW; i++) fprintf(pFile, "%010d\r\n", ROWIND[i]);
		//block of data indices
		for (int i = 0; i<NNZERO; i++) fprintf(pFile, "%020.16f\r\n", VALUES[i]);

		fclose(pFile);
	}
	catch (const exception e) {
		cout << "Exception saving HarwellBoeing file: ";
		return -1;
	}
	

	return 0;
}



int HarwellBoeing::Load(string filename, double Matrix[]){ return 0; }
