#ifndef _ANALYSIS_H
#define _ANALYSIS_H

#include "../mode_null/null_data.h"

#include <boost/math/distributions/hypergeometric.hpp>

using namespace std;
using namespace boost::math;


#define RESAMPLING 0
#define NO_RESAMPLING 1 

struct genomic_region{
    std::string chrom;
    int start, end;
};

struct positions{
    int start, end;
};

class analysis_cpp {
public:

    // ANALYSIS VARIABLES // 
    int window_size;
    float window_maf;
    int random_variants;
    unsigned int mode;

    // QTL 
    int qtl_count;
    vector < string >  qtl_id;
    vector < int > qtl_dist_phe_var;
    vector < string > qtl_phe_strand;
    vector < string > qtl_chr;
    vector < int > qtl_start;
    vector < int > qtl_end;
    vector < float > qtl_maf;

    unordered_map < string, int > qtl_map;


    // NULL
    int null_count;
    vector < string > null_id;
    vector < string > null_chr;
    vector < int > null_start;
    vector < int > null_end;
    vector < float > null_maf;
    vector < int > upstream_distance;
    vector < int > downstream_distance; 
    unordered_map <string, float> map_maf;
    vector < unsigned int > nominal;

    // NULL DISTRIBUTION
    vector < string > nulldistribution; // Contains all the variants for the null distribution for the functional enrichment analysis
    vector < genomic_region > nulldistribution_regions;
    unordered_map <string, int> removeVar; 
    int empty=0;
    int below_random_threshold =0;
    // NULL DISTRIBUTION CREATION VARIABLES 
    unsigned int seed=12345;

    // READ INTERSECTION DATA
    int intersection_count;
    set < string > peaks;
    vector < string > intersection_snp;
    vector < string > intersection_peaks;

    // READ PHENOTYPES
    int phen_count;
    std::vector < string > phen_chr;
    std::vector < int > phen_start;
    std::vector < int > phen_end;
    
    unordered_map < std::string, positions > phen_index;

    // READ DATA 
    void readQTL(string);
    void readNull(string);
    void readIntersection(string);
    void readPhenotypes(std::string);

    // PROCESS
    void createNullDistribution(string);
    void functionalEnrichment(string);
    //void filterQTL();
    void performIntersect(string);
    bool removeVarFromNull(int&);

    // INLINE FUNCTIONS FOR PERFORMANCE
    double fisher_test(unsigned, unsigned, unsigned, unsigned);
    double odds_ratio(unsigned, unsigned, unsigned, unsigned);
    int findSNP(vector<string>,string);
    bool has_chr(std::string);

    // OPTIONS 
    boost::program_options::options_description option_descriptions;
	boost::program_options::variables_map options;
};

//***************************************************************//
//******************** DECLARE FUNCTIONS *************************//
//***************************************************************//
void analysis_main(vector < string > & argv);



//***************************************************************//
//******************** INLINE FUNCTIONS *************************//
//***************************************************************//

inline bool analysis_cpp::has_chr(std::string chr)
{
    const char* chrom = "chr";
    bool b = boost::algorithm::contains(chr, chrom);
    return b;
}


inline int analysis_cpp::findSNP(vector < string > vec,string snp){

    vector<string>::iterator itr = find(vec.begin(), vec.end(), snp);
    if(itr != vec.end()){
        return (distance(vec.begin(),itr));
    }else{
        cout << "Element not found in vector" << endl;
        return(-1);
    }
}

inline double analysis_cpp::fisher_test(unsigned a, unsigned b, unsigned c, unsigned d) {
		unsigned N = a + b + c + d;
		unsigned r = a + c;
		unsigned n = c + d;
		unsigned max_for_k = min(r, n);	
		unsigned min_for_k = (unsigned)max(0, int(r + n - N));
		hypergeometric_distribution<> hgd(r, n, N);	
		double cutoff = pdf(hgd, c);
		double tmp_p = 0.0;
		for(int k = min_for_k;k < max_for_k + 1;k++) {
				double p = pdf(hgd, k);
				if(p <= cutoff) tmp_p += p;
		}
		return tmp_p;
}

inline double analysis_cpp::odds_ratio(unsigned a, unsigned b, unsigned c, unsigned d){
    return(((double)a*(double)d) / ((double)b*(double)c));
}


#endif /* _FENRICH_CPP_H */
