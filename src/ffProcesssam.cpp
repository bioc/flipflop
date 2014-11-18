#include <iostream>
#include <cstring>
#include "samio.h"
#include "common.h"

extern "C" {  
  void ffProcesssam(char ** in, char ** pre, char ** an, char ** pair, char ** minr, char ** minc, char ** minj, char ** verb){
    string inSamFile = string(in[0]);  
    string prefix = string(pre[0]);
    string paired = string(pair[0]);
    string verbose = string(verb[0]);
    string minReadNum = string(minr[0]);
    string minJuncCount = string(minj[0]);
    string minCvgCut = string(minc[0]);
    vector<string> cline;  
    /* processsam activation */
    cline.push_back(string("processsam"));
    /* consider reads as single-end */
    if(paired=="FALSE"){
       cline.push_back(string("--single-only"));
    }
    /* prefix name */
    cline.push_back(string("--prefix"));
    cline.push_back(prefix);
    /* verbosity */
    cline.push_back(string("--verbose"));
    cline.push_back(verbose);
    /* optional annotation */
    if(strcmp("", an[0])){
      //cline.push_back(string("--annotation"));
      cline.push_back(string("-x"));
      cline.push_back(string(an[0]));
    }
    /* minimum number of reads in a gene */
    cline.push_back(string("--min-read-num"));
    cline.push_back(minReadNum);
    /* minimum number of reads for considering a valid junction */
    cline.push_back(string("--min-junc-count"));
    cline.push_back(minJuncCount);
    /* coverage cut off for segmenting sub-exons */
    cline.push_back(string("--min-cvg-cut"));
    cline.push_back(minCvgCut);
    cline.push_back(inSamFile);
    /* initialize the boring global variables */
    initGlobalVariables();
    /* start the job */
    readSamFile(inSamFile, prefix, cline);  
  }
}
