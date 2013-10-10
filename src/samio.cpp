#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <algorithm>
#include <climits>

#include "samio.h"
#include "parseopt.h"
#include "auxiliaryio.h"

//#include "readgroup.h"

inline long abs2(long a){
  return a>0?a:-a;
}

/* 
Initialize the boring and tiring global variables.
These variables are the one in common and structdef,
but the value are changed when changing the option (annotation for instance). 
Hence we forced it to have the default value at the beginnin of ffProcesssam, 
before to read the options given by the user. 
*/
void initGlobalVariables(){
    n_INST = 0;
    REFONLY=false;
    FIXRANGE=false;
    FIXBOUND=false;
    VERBOSE=0;
    MIN_CVG_FRACTION=0.25;
    MINEXONCVG=2;
    OUTPUT_INDIVIDUAL_COVERAGE=true;
    DEFAULT_MIN_JUNCTION=1;
    MIN_GRANGE_DISTANCE=100;
    MIN_GRANGE_READ_CNT=4;
    MAX_INSTANCE=-1;
    MAX_PE_DISTANCE=700000;
    MAX_N_SEGS=500000;
    USE_SINGLE_ONLY=0;
    MAX_READ_SPAN=50000;
    MAX_READ_OVER=100000;
    //MIN_SEG_FL_OVERLAP=10;
    MIN_SEG_FL_OVERLAP=0;
    STRANDED_RNASEQ=0;
    /* provided gene ranges. each range may contain multiple segments */
    C_GENE_RANGE.clear();
    /* Provided junction boundaries */
    C_JUNCTION.clear();
    /*Provided annotation*/ 
    C_ANNO.resetAll();
    C_RANGE_ANNO.resetAll();
}


/*
analyze sam file.
Return value: 0 if all reads are single-ended, 1 if some are paired-ended.
*/
int readSamFile(string inSamFile, //input file
    vector<string> args // we need to print the command line into the instance file
    ){

  //Parsing arguments and set up parameters
  if (parseopt(args)==-1){
    return -1;
  }

  vector<string> rangechrs;
  vector<range_t > rangerange;
  bool outinstance=true;
  //Preparing output files
  prepareAuxFile(args, inSamFile, outinstance);

  ifstream fifs;
  if(inSamFile!="-"){
    fifs.open(inSamFile.c_str());
    if(!fifs.is_open()){
      // cerr<<"Error opening input file "<<inSamFile<<endl;
      return -1;
    }
  }
  istream &ifs= (inSamFile=="-")?cin:fifs;
  // cerr<<"Input file: "<<inSamFile<<endl;

  int appearpereads=-1;
  if(USE_SINGLE_ONLY==1)appearpereads=0;
  //saving flag and cigar for pair-end 
  map<range_t, int> flagmap;
  map<range_t, string> cigarmap;
  map<range_t, int> paircount;

  string prevrname="";
  range_t dp=make_range_t(0,0);
  range_t dp2=dp;
  
  range_t currentrange=dp;
  long currentrcount=0;
  int rangecounter=1;


  int instanceid=0; //instance id
  int currentreadlen=0;
  long linecount=0;
  int totalnumread=0; // total number of reads in the SAM file 
  int totalpair=0; // total number of pairs in the SAM file 

  Align calign;
  ReadGroup rd;

  //----------------The main loop of parsing SAM file-------------------------
  while(true){

    string oneline;
    getline(ifs,oneline);
    linecount++;
    if(ifs.eof())break;
    if(oneline[0]=='@') continue;
    if(calign.parse(oneline)==-1){
      // cerr<<"Error parsing SAM records at line "<<linecount<<endl;
    }
    bool usesingle=USE_SINGLE_ONLY;
    //if the paired-end distance is too large, treat them as single-end reads
    if(calign.isPairedEnd()) {
       totalpair++;  // ATTENTION ELSA 17 SEPT 13
       if( abs2(calign.plen)>MAX_PE_DISTANCE)
         usesingle=true;
       //do not support mapping to different chromosomes now
       if(calign.rnext!="*" && calign.rnext!="=")
         usesingle=true;
    }
    if(!calign.isValid())
      continue;
  
    //write to read info
    writeoneline2readinfo(calign);
    if(linecount==1){
      ReadGroup::setStatReadLen(calign.getReadLen());
      ReadGroup::setStatChr(calign.rname);
    }
    if(oneline[0]!='@') totalnumread++;
    //write to instance, if possible
    //update and write generange file
    if(calign.rname!=prevrname || calign.pos>currentrange.second+MIN_GRANGE_DISTANCE){
  
      rd.clearPairInfo();   
      if(currentrange.second!=0)
        write2rangeandinstance(rd,currentrange);
      rd.clear();
      //If gene ranges are fixed, get the next gene range
      if(FIXRANGE){
        string crname=calign.rname;
        int getnextret=0;
        
        while((getnextret=C_GENE_RANGE.getNext(calign.rname,calign.pos,crname,currentrange))==2){
          // if annotation's gene range < read's gene range, output annotation gene range
          ReadGroup::setStatChr(crname);
          rd.setRange(currentrange);
          write2rangeandinstance(rd,currentrange);
          // cout<<"\nJumping range, read: "<<calign.rname<<", range: "
          //  <<calign.rname<<":"<<currentrange.first<<"-"<<currentrange.second<<endl;
        }
        if(getnextret==-1)
           break;
        ReadGroup::setStatChr(crname);
        rd.setRange(currentrange);
        // if(currentrange.first>0)
          // cout<<"\nFix range, read: "<<calign.rname<<", range: "
          //  <<calign.rname<<":"<<currentrange.first<<"-"<<currentrange.second<<endl;
      }else{
        //update current range and current r count
        currentrange=calign.getRange(usesingle);
      }
      //update current range w.r.t. existing ref isoforms
      //this works no matter whether FIXRANGE is true
      //C_ANNO.checkOverlapRange(currentrange,rname);
    }//end if(pos>currentrange.second+MIN_GRANGE_DISTANCE)
      
    //update the current range
    if(FIXRANGE){
      range_t rrt=calign.getRange(usesingle);
      if(rrt.first>currentrange.first && rrt.second<currentrange.second)
        rd.add(calign);
    }else{
      if(calign.pos<currentrange.first){
          // cerr<<"Error: the SAM file MUST be sorted!\n";
      }
      else{
        currentrcount++;
        long antd=calign.getRange(usesingle).second;
        if(antd>currentrange.second){
          currentrange.second=antd;
        }
        //this works no matter whether FIXRANGE is true
        C_ANNO.checkOverlapRange(currentrange,calign.rname);
      }
      rd.add(calign);
    }
    prevrname=calign.rname;
    if(calign.rname!=prevrname)
      ReadGroup::setStatChr(calign.rname);
    ReadGroup::setStatChr(calign.rname);
  
    if(MAX_INSTANCE!=-1 && instanceid>MAX_INSTANCE)
      break;
  
    if(rd.size()>10000 && rd.size()%10000==1){
      // cout<<"--Reads: "<<rd.size()<<", range:["<<currentrange.first<<","<<currentrange.second<<"]"
      //    <<", current read: "<<calign.pos<<endl;
    }
  }//end of iteration of ifs
  //write the last one
  write2rangeandinstance(rd,currentrange);
  //write empty ref seq ranges, if possible
  currentrange.first=LONG_MAX;
  rd.clear();  
  if(FIXRANGE){
    // cout<<"\nEnd of read input, writing remaining instances...\n";
    while(!C_GENE_RANGE.isEnd()){
      currentrange=C_GENE_RANGE.getRange();
      ReadGroup::setStatChr(C_GENE_RANGE.getChr());
      rd.setRange(currentrange);
      write2rangeandinstance(rd,currentrange);
      C_GENE_RANGE.inc();
    }
  }
  //write the total number of reads in a separate file
  ofstream fichier("totalnumread", ios::out | ios::trunc);
  fichier<<"@Total Number of Reads\n"<<totalnumread<<"\n";
  fichier<<"@Total Number of Pairs\n"<<totalpair/2;
  fichier.close();
  fifs.close();
  closeAuxFile();

  return appearpereads;
}


