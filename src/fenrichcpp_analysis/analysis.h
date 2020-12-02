#ifndef _ANALYSIS_H
#define _ANALYSIS_H

#include "../fenrichcpp.h"

using namespace std;

class analysis_cpp {
public:

    // ANALYSIS VARIABLES // 
    int window_size;
    float window_maf;
    int random_variants;

    // QTL 
    int qtl_count;
    vector < string > qtl_id;
    vector < int > qtl_dist_phe_var;
    vector < string > qtl_phe_strand;
    vector < int > qtl_dist_phe_var_from;
    vector < int > qtl_dist_phe_var_to;
    vector < float > qtl_maf;
    vector < float > qtl_maf_from;
    vector < float > qtl_maf_to;

    // NULL

    int null_count;
    vector < string > null_id;
    vector < float > maf;
    vector < int > upstream_distance;
    vector < int > downstream_distance; 


    // READ DATA 
    void readQTL(string);
    void readNull(string);
    



}
