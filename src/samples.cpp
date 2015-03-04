#include <iostream>
#include <sstream>
#include <cstdio>
#include <fstream>
//#include <assert.h>
#include "samples.h"

#include <R.h> // add that for using Rprintf

const Samples* Samples::instance;

// in place trimming
void trim(std::string& str)
{
   size_t pos = str.find_first_not_of(" \t");
   if (pos != std::string::npos) {
     str.erase(0, pos);
   }
   pos = str.find_last_not_of(" \t");
   if (pos != std::string::npos) {
      str.erase(pos+1);
   }
}


Samples::Samples(){
 //assert(!instance);
 instance=this; 
}
int Samples::readsamples(const string& filename) // Samples constructor
{
   ifstream fifs;
   fifs.open(filename.c_str());
   if(!fifs.is_open()){
     Rprintf("Error opening input sample file %s\n",filename.c_str());
     return -1;
   }
   Rprintf("Input sample file %s\n",filename.c_str());

   while(!fifs.eof()){
    string oneline;
    getline(fifs,oneline);  
    trim(oneline); // trim the sample name in case there is a blank
    if(oneline.length()>0){
      samples.push_back(oneline);
    }
   }
   return 0;
}

#if 0 
int main(int argc, char* argv[]){
  Samples s1; 
  for(int n=1;n<argc;n++){
    if(s1.readsamples(argv[n])<0){
      return 1;
    }
  }
  return 0;
}
#endif
