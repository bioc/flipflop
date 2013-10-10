/**
SAMIO.H
**/

#ifndef SAMIO_H
#define SAMIO_H

#include <string>
#include <vector>
#include "structdef.h"

using namespace std;

//int ffProcesssam(char * in, char * pre, char * an);

/*
Main entry for analyzing sam file
read sam file, and generate read mapping file for isoinfer
*/

void initGlobalVariables();

int readSamFile(string inSamFile, 
		vector<string> args
		);


#endif
