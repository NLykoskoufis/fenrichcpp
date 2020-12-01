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
    vector < string > qtl_chr;
    vector < int > qtl_start;
    vector < int > qtl_dist_phe_var;
    vector < string > qtl_phe_strand;


    // NULL

    int null_count;
    vector < string > null_id;
    vector < float > maf_from;
    vector < float > maf_to;
    vector < int > upstream_from; 
    vector < int > upstream_to;
    vector < int > downstream_from;
    vector < int > downstream_to;

    



}
