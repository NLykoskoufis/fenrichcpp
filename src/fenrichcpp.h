
#ifndef _FENRICH_CPP_H
#define _FENRICH_CPP_H


#include "compression_io.h"

//INCLUDE STANDARD TEMPLATE LIBRARY USEFULL STUFFS (STL)
#include <vector>
#include <numeric>
#include <string>
#include <cmath>

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

using namespace std;

class fenrich_cpp {
public:
    // REGIONS 
    string regionPhenotype_chr = "NA";
    string regionPhenotype;

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
    int genotype_count=0;
    vector < string > genotype_id;
    vector < string > genotype_chr;
    vector < int > genotype_start;
    vector < int > genotype_end;
    vector < vector < float > > genotype_val;
    
    // NOMINAL QTL 
    int qtl_count;
    vector < string > qtl_id;

    // MANAGEMENT 
    void fenrichcpp_createTEnull(string);

    // CONSTRUCTOR / DESTRUCTOR 
    //fenrich_cpp();
    //~fenrich_cpp();
    //void clear();

    // READ DATA 
    void readPhenotypes(string);
    void readGenotypes(string);
    void readSignificantQTL(string);

    // COMPUTATION METHODS [ALL INLINE FOR SPEED]
    float getMedian(vector<float>);
    float getMean(vector <float>);
    float getMAF(vector <float>);
};


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





