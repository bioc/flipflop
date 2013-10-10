#include "parseopt.h"
#include "auxiliaryio.h"
#include <cstdlib>
#include <iostream>


/*
Parsing the parameters, setting up variables and structures
*/
int parseopt(vector<string> args){

  // analyze input arguments
  for(int i=0;i<args.size();i++){
    //these are parameters that need to process first
    if(string(args[i])=="-b" || string(args[i])=="--single-only"){
      USE_SINGLE_ONLY=1;
      // cout<<"Use single-end reads only.\n";
    }
    //if(string(args[i])=="-r" ||string(args[i])=="--ref-only" ){
      //depreciated. Using -x option will automatically turn REFONLY on.
      //REFONLY=true;
    //}
    if(string(args[i])=="--verbose"){
       VERBOSE=atoi(args[i+1].c_str());
    }   
    if(string(args[i])=="-v" || string(args[i])=="--no-coverage"){
      OUTPUT_INDIVIDUAL_COVERAGE=false;
    }
    if(i<args.size()-1){
      if(string(args[i])=="-g" || string(args[i])=="--min-gap-length"){
        MIN_GRANGE_DISTANCE=atoi(args[i+1].c_str());
      }
      else if(string(args[i])=="-c" || string(args[i])=="--min-read-num"){
        MIN_GRANGE_READ_CNT=atoi(args[i+1].c_str());
      }
      else if(string(args[i])=="-k" || string(args[i])=="--max-pe-span"){
        MAX_PE_DISTANCE=atoi(args[i+1].c_str());
        // cout<<"Maximum pair-end span: "<<MAX_PE_DISTANCE<<endl;
      }
      else if(string(args[i])=="-s" || string(args[i])=="--max-num-instance"){
        MAX_INSTANCE=atoi(args[i+1].c_str());
      }
      else if(string(args[i])=="-e"|| string(args[i])=="--segment-bound"){
        //MAX_N_SEGS=atoi(args[i+1].c_str());
        //cout<<"Maximum number of segments: "<<MAX_N_SEGS<<endl;
        FIXBOUND=true;
        string juncfile=string(args[i+1]);
        // cout<<"Junction file: "<<juncfile<<endl;
        C_JUNCTION.clear();
        readgenerangefile(juncfile,C_JUNCTION);
      }
      else if(string(args[i])=="-r" || string(args[i])=="--range"){
        //read gene range file
        string annobedfile=string(args[i+1]);
        //providing existing gene annotation
        if(annobedfile!=""){
          if(parseannobed(annobedfile,C_RANGE_ANNO)==-1) return -1;

          //boundary is fixed
          FIXRANGE=true;
          C_GENE_RANGE.clear();
          vector<string> allchr=C_RANGE_ANNO.getChrom();
          for(int ai=0;ai<allchr.size();ai++){
            vector<range_t> allr=C_RANGE_ANNO.getCluster(allchr[ai]);
            for(int ri=0;ri<allr.size();ri++){
              //add gene range
              C_GENE_RANGE.push_back(allchr[ai],allr[ri]);
            }
          }
        }
 
      }
      else if(string(args[i])=="-j" || string(args[i])=="--min-junc-count"){
        DEFAULT_MIN_JUNCTION=atoi(args[i+1].c_str());
        // cout<<"Minimum junction count: "<<DEFAULT_MIN_JUNCTION<<endl;
      }
      else if(string(args[i])=="-x" || string(args[i])=="--annotation"){
        string annobedfile=string(args[i+1]);
        // cout<<"Bed file: "<<annobedfile<<endl;
        //providing existing gene annotation
        if(annobedfile!=""){
          if(parseannobed(annobedfile,C_ANNO)==-1) return -1;

          //boundary is fixed
          FIXBOUND=true;
          //currently, only use annotation to build profiles
          REFONLY=true;
          if(REFONLY){
            //if this option is on, fill the FIXRANGE variables.
            //this will automatically turn the variable FIXRANGE=true in later codes
            //also, fill the C_JUNCTION
            C_GENE_RANGE.clear();
            vector<string> allchr=C_ANNO.getChrom();
            for(int ai=0;ai<allchr.size();ai++){
              vector<range_t> allr=C_ANNO.getCluster(allchr[ai]);
              for(int ri=0;ri<allr.size();ri++){
                //add gene range
                C_GENE_RANGE.push_back(allchr[ai],allr[ri]);
                //add exon-intron boundaries
                ReadGroup rg=C_ANNO.getReadGroup(allchr[ai],allr[ri]);
                vpos_t &cs=rg.s(); vpos_t& ce=rg.e();
                for(int i=0;i<cs.size();i++){
                  //are they sorted?
                  for(int j=0;j<cs[i].size();j++){
                    long cl=cs[i][j], cr=ce[i][j];
                    if(cl!=allr[ri].first)
                      C_JUNCTION.push_back(allchr[ai],make_range_t(cl,cl));
                    if(cr!=allr[ri].second)
                      C_JUNCTION.push_back(allchr[ai],make_range_t(cr+1,cr+1));//the real junction position should be 1 base away to the end position
                  }
                }
              }
            }
          }//end if REFONLY
        }
      }
      else if(string(args[i])=="-u" || string(args[i])=="--min-cvg-cut"){
        MIN_CVG_FRACTION=atof(args[i+1].c_str());
        // cout<<"Minimum coverage cutoff: "<<MIN_CVG_FRACTION<<endl;
      }
      else if(string(args[i])=="-f"){//depreciated
        //string boundfile=string(args[i+1]);
        //cout<<"Gene range file: "<<boundfile<<endl;
        //C_GENE_RANGE.clear();
        //readgenerangefile(boundfile,C_GENE_RANGE);
      }
      else if(string(args[i])=="-d" || string(args[i])=="--direction"){
        string nsr=string(args[i+1]);
        if(nsr=="+")
          STRANDED_RNASEQ=1;
        if(nsr=="-")
          STRANDED_RNASEQ=-1;
      }
    }
  }
  
  //cout<<"Minimum gap length for two genes (mingapdist): "<<MIN_GRANGE_DISTANCE<<endl;
  //cout<<"Minimum number of reads in one gene (mingrangecount): "<<MIN_GRANGE_READ_CNT<<endl;
 
  if(C_GENE_RANGE.size()>0){
    FIXRANGE=true;
    MIN_GRANGE_DISTANCE=0;
    MIN_GRANGE_READ_CNT=0;
  }

  // cout<<"Minimum gap length for two genes (mingapdist): "<<MIN_GRANGE_DISTANCE<<endl;
  // cout<<"Minimum number of reads in one gene (mingrangecount): "<<MIN_GRANGE_READ_CNT<<endl;
  
  // if(VERBOSE==1){
    // cout<<"Using fix range: "<<FIXRANGE<<",size:"<<C_GENE_RANGE.size()<<endl;
  // }
  C_JUNCTION.sort();
  
  if(C_JUNCTION.check()==false){
    // cout<<"Check failed.\n";
    return -1;
  }
  return 0;


}


