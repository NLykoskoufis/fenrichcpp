
#ifndef _FENRICH_CPP_H
#define _FENRICH_CPP_H


#include "../../lib/compression_io.h"
#include "../../lib/Instrumentor.h"

//INCLUDE STANDARD TEMPLATE LIBRARY USEFULL STUFFS (STL)
#include <vector>
#include <numeric>
#include <string>
#include <cmath>
#include <map> 
#include <iterator>
#include <unordered_map>
#include <limits>
#include <set> 
#include <algorithm>
#include <random> 
#include <iomanip>

//INCLUDE BOOST LIBRARIES
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

#include <htslib/hts.h>
#include <htslib/kseq.h>
#include <htslib/sam.h>

extern "C" {
	#include <htslib/vcf_sweep.h>
	#include <htslib/synced_bcf_reader.h>
	#include <htslib/vcf.h>
	#include <htslib/vcfutils.h>
}

#define PROFILING 1
#if PROFILING
#define PROFILE_SCOPE(name) InstrumentationTimer timer##__LINE__(name)
#define PROFILE_FUNCTION() PROFILE_SCOPE(__FUNCSIG__)
#else
#define PROFILE_SCOPE(name)
#endif


using namespace std;

class fenrich_cpp {
public:
    
    string genomic_region;

    // PHENOTYPES 
    int phenotype_count;
    vector < string > phenotype_id;
    vector < string > phenotype_chr;
    vector < int > phenotype_start;
    vector < int > phenotype_end;
    vector < string > phenotype_strand;
    vector < vector < float > > phenotype_val;
    vector < string > phenotype_length;
    vector < string > phenotype_neg;
    
    // GENOTYPES 
    int genotype_count;
    vector < string > genotype_id;
    vector < string > genotype_chr;
    vector < int > genotype_start;
    vector < int > genotype_end;
    vector < vector < float > > genotype_val;
    

    // NOMINAL QTL 
    int qtl_count;
    unordered_map<string,unsigned int> nomQTL;

    //DATA REGION
	bool setPhenotypeRegion(string);
	bool setGenotypeRegion(string);
	void setPhenotypeRegion(int, int);

    // MANAGEMENT 
    void fenrichcpp_createTEnull(string);

    // READ DATA 
    void readPhenotypes(string);
    void readGenotypes(string);
    void readSignificantQTL(string);

    // COMPUTATION METHODS [ALL INLINE FOR SPEED]
    float getMedian(vector<float>);
    float getMean(vector <float>);
    float getMAF(vector <float>);

    // OPTIONS 
    boost::program_options::options_description option_descriptions;
	boost::program_options::variables_map options;

};

//***************************************************************//
//******************** DECLARE FUNCTIONS *************************//
//***************************************************************//
void null_main(vector < string > &);


//***************************************************************//
//******************** INLINE FUNCTIONS *************************//
//***************************************************************//

inline float fenrich_cpp::getMAF(vector < float > gt){
  float size = gt.size();
  if (size ==0){
    return -1;
  }
  else{
    float AF=0;
    for(int i = 0; i < size; i++){
      if(gt[i] == 0){
        AF = AF + 2;
      }else if(gt[i] == 1){
        AF = AF +1;
      }else{
        continue;
      }
    }
    if((AF / (2*size)) > 0.5 ){
      return (1-(AF/(2*size)));}
    else{
      return (AF/(2*size));
    }
  }
}


inline float fenrich_cpp::getMedian(vector<float> scores)
{
  size_t size = scores.size();
  if (size == 0)
  {
    return 0;  // Undefined, really.
  }
  else
  {
    sort(scores.begin(), scores.end());
    if (size % 2 == 0)
    {
      return (scores[size / 2 + 1] + scores[size / 2]) / 2;
    }
    else 
    {
      return scores[round(size / 2)];
    }
  }
}

inline float fenrich_cpp::getMean(vector <float> scores){
    float sum = std::accumulate(scores.begin(),scores.end(),0.0);
    return sum / scores.size();
}



#endif /* _FENRICH_CPP_H */





