/*
 * read sample names from file
*/

#include <iostream>
#include <vector>
#include <string>
using namespace std;

// 2015-01-15 ELSA
class Samples
{
  vector<string> samples;
  static const Samples* instance;
  public:
    Samples();
    int readsamples(const string& filename);
    const vector<string>& getSamples() const{
      return samples;
    } 
    static const Samples* getInstance(){return instance;}
};

// in place triming
extern void trim(std::string& str);
