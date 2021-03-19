#include "analysis_data.h"


void analysis_cpp::performIntersect(std::string fout){
    //output_file fdo("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/eqtls/enrichment/analysis_functionalEnrichment/fenrich_cpp/test/fenrich_cpp_overlapped.txt");

    
    
    int qtl_overlap=0;
    int null_overlap=0;
    for(int i=0; i<qtl_id.size(); i++)
    {   
        std::cout << qtl_chr[i] << std::endl;
        unordered_map<std::string,positions>::iterator got = phen_index.find(qtl_chr[i]);
        if (got == phen_index.end()) continue;
        for(int j=got->second.start; j < got->second.end; j++)
        {
            //std::cout << i << std::endl;
            if(qtl_chr[i].compare(phen_chr[j]) != 0) std::cout << "Problem" << std::endl;
            if(qtl_start[i] <= phen_end[j] && qtl_end[i] >= phen_start[j])
            {
                //fdo << qtl_id[i] << std::endl;
                qtl_overlap++;
                break; // Count only once!
            }
        }
    }
    
    
    for(int i=0; i < nulldistribution.size();i++)
    {   
        unordered_map<std::string,positions>::iterator got = phen_index.find(nulldistribution_regions[i].chrom);
        if (got == phen_index.end()) continue;
        for(int j=got->second.start; j < got->second.end; j++)
        {
            if(nulldistribution_regions[i].chrom.compare(phen_chr[j]) != 0) std::cout << "Problem" << std::endl;
            if(nulldistribution_regions[i].start <= phen_end[j] && nulldistribution_regions[i].end >= phen_start[j])
            {
                null_overlap++;
                break;
            }
        }
    }
    int qtl_no_overlap = qtl_id.size() - qtl_overlap;
    int null_no_overlap = nulldistribution.size() - null_overlap;
    
    if(qtl_overlap == 0 || null_overlap == 0)
    {
        qtl_overlap++;
        null_overlap++;
    } 

    std::cout << qtl_overlap << " " << qtl_no_overlap <<  std::endl;
    std::cout << null_overlap << " " << null_no_overlap << std::endl;

    double fisher = fisher_test(qtl_overlap,null_overlap, qtl_no_overlap,null_no_overlap);
    double oddsratio = odds_ratio(qtl_overlap,null_overlap, qtl_no_overlap,null_no_overlap);

    std::cout << " ** Performing Enrichment using Fisher exact test." << std::endl;
    
    output_file fdo (fout);
    fdo << "A\tB\tC\tD\tOddsRatio\tpvalue" << endl;
    
    cout << "qtl overlap: " << qtl_overlap << "\nqtl no overlap: " << qtl_no_overlap << "\nnull_overlap: " << null_overlap << "\nnull_no_overlap: " << null_no_overlap << endl;
    cout << fisher << endl;
    cout << oddsratio << endl;
    std::ostringstream strs;
    strs << fisher;


    fdo << to_string(qtl_overlap) << "\t" << to_string(null_overlap) << "\t" << to_string(qtl_no_overlap) << "\t" << to_string(null_no_overlap) << "\t" << to_string(oddsratio) << "\t" << strs.str() << endl;

}