#include "analysis_data.h"

/* Here I am just reading Halit's chipseq_anno.txt and trying to see whether his and mine are giving different results""
"
This is an alternative to overlapping the peaks with the eQTLs. Here we just need to check whether the SNPs are overlapping with the TFs.""
*/

void analysis_cpp::readPeakAnnotation(std::string fanno){

    input_file fd (fanno);
    std::string buffer; 
    std::vector < std::string > line;

    int linecount = 0;
    while(getline(fd, buffer))
    {   
        if (linecount % 1000000 == 0) cout << std::to_string(linecount) << " lines read" << endl;
        linecount++;
        boost::split(line, buffer, boost::is_any_of("\t"));
        std::vector < std::string > anno; 
        
        boost::split(anno, line[1], boost::is_any_of(";"));
        if(anno.size() == 1) continue;
        std::set < std::string > anno_set;
        for(int p=1; p < anno.size(); p++)
        {   
            anno_set.insert(anno[p].substr(15,anno[p].size()));
        }
        anno_dico.insert(std::make_pair(line[0],anno_set));
    }
    fd.close();
    std::cout << "Read " << anno_dico.size() << " overlapping with phenotype data" << std::endl;
}


void analysis_cpp::performEnrichment(std::string fout){
    std::cout << "TODO" << std::endl;
    std::cout << "Performing Enrichment for " << mark << std::endl;

    int qtl_overlap =0, null_overlap = 0;
    for(int i=0; i<qtl_id.size(); i++)
    {   
        std::unordered_map < std::string, std::set < std::string > >::iterator got = anno_dico.find(qtl_id[i]);
        if (got == anno_dico.end()) continue;
        if (got->second.find(mark) != got->second.end()) qtl_overlap++;
    }

    for(int i=0; i<nulldistribution.size(); i++)
    {   
        std::unordered_map < std::string, std::set < std::string > >::iterator got = anno_dico.find(nulldistribution[i]);
        if (got == anno_dico.end()) continue;
        if (got->second.find(mark) != got->second.end()) null_overlap++;
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


