#include "analysis_data.h"

using namespace std;

void analysis_cpp::createNullDistribution(string fout){


    PROFILE_FUNCTION();

    /*
    For each eQTL:
        1. Loop over all null and get all the SNPs matching MAF, distance and that are not eQTLs (nominal != 1).
        2. From vector randomly select 10 of them.
        3. For the next eQTL, check that the eQTLs included in vector are not the same as the previous one. This step might be long maybe its better to save these ones inside a dictionary for faster search than a vector. 

    When the null distribution of SNPs are done, then we can perform the functional enrichment. 
    */
   std::string output_name = fout + "_null_distribution_used";
   output_file fdo (output_name);
   fdo << "qtl_id\tmaf\tdist_phe_var\tmaf_from\tmaf_to\tdist_from\tdist_to\trandom_variants" << endl;
   for(int i=0; i<qtl_id.size(); i++){
        vector < string > toRandomPeak;
        vector < genomic_region > toRandomPeak_regions;
        cout << "Processed " << i+1 << " eQTLs." << endl;
        
        // GET UPSTREAM AND DOWNSTREAM 
            // If d > 0 --> Downstream and if from is smaller than 0 meaning that bin is tresspassing into upstream territory, set "to" to 0. We do not want to go upstream if variant is downstream. The same for second if but for downstream.
        float qtl_maf_from;
        float qtl_maf_to; 
        int qtl_dist_phe_var_from;
        int qtl_dist_phe_var_to;

        qtl_maf_from = qtl_maf[i] - (qtl_maf[i] * window_maf);
        qtl_maf_to = qtl_maf[i] + (qtl_maf[i] * window_maf);

        if(qtl_dist_phe_var[i] < 0){
            qtl_dist_phe_var_to = (qtl_dist_phe_var[i] + window_size) > 0 ? 0:(qtl_dist_phe_var[i] + window_size);
            qtl_dist_phe_var_from = qtl_dist_phe_var[i] - window_size;
        }else{
            qtl_dist_phe_var_from = (qtl_dist_phe_var[i] - window_size) < 0 ? 0:(qtl_dist_phe_var[i] - window_size);
            qtl_dist_phe_var_to = qtl_dist_phe_var[i] + window_size;
        }
        // Writing file containing distances and randomly chose variants.
        fdo << qtl_id[i] << "\t" << to_string(qtl_maf[i]) << "\t" << to_string(qtl_dist_phe_var[i]) << "\t" << to_string(qtl_maf_from) << "\t" << to_string(qtl_maf_to) << "\t" << to_string(qtl_dist_phe_var_from) << "\t" << to_string(qtl_dist_phe_var_to) << "\t";

        for(int s=0; s < null_count; s++){
           // Check whether eQTL dist and MAF windows are matching with any TEs. 
           // If they do not match, through message
           
           if(qtl_dist_phe_var[i] > 0){
               if(qtl_dist_phe_var_from <= downstream_distance[s] && qtl_dist_phe_var_to>= downstream_distance[s]){
                   if(qtl_maf_from <= null_maf[s] && qtl_maf_to >= null_maf[s]){
                       if(nominal[s] != 1 && find(nulldistribution.begin(),nulldistribution.end(), null_id[s]) == nulldistribution.end()){
                           toRandomPeak.push_back(null_id[s]);
                           toRandomPeak_regions.push_back({null_chr[s], null_start[s], null_end[s]});
                       }
                   } 
               }
           }else{
               if(qtl_dist_phe_var_from <= upstream_distance[s] && qtl_dist_phe_var_to>= upstream_distance[s]){
                    if(qtl_maf_from <= null_maf[s] && qtl_maf_to >= null_maf[s]){
                        if(nominal[s] != 1 && find(nulldistribution.begin(),nulldistribution.end(),null_id[s]) == nulldistribution.end()){
                            toRandomPeak.push_back(null_id[s]);
                            toRandomPeak_regions.push_back({null_chr[s], null_start[s], null_end[s]});
                        }
                    }
                }
            }
        }// NULL for loop
        if(toRandomPeak.size() == 0){
            cout << "No variants found" << endl;
            empty++;

        }else{
            //cout << toRandomPeak.size() << endl;
            std::vector < int > v(toRandomPeak.size());
            std::iota(std::begin(v), std::end(v), 0);
            std::shuffle(v.begin(), v.end(), std::default_random_engine(seed));
            
            //std::shuffle(toRandomPeak.begin(), toRandomPeak.end(), std::default_random_engine(seed));
            if(toRandomPeak.size() >= 10){
                
                for(int r=0; r < 10; r++){
                    nulldistribution.push_back(toRandomPeak[v[r]]);
                    nulldistribution_regions.push_back(toRandomPeak_regions[v[r]]);
                    fdo << toRandomPeak[v[r]];
                    if(r < 10) fdo << ",";
                }
                fdo << endl;
            }else{
                
                for(int r=0; r < toRandomPeak.size(); r++){
                    nulldistribution.push_back(toRandomPeak[v[r]]);
                    nulldistribution_regions.push_back(toRandomPeak_regions[v[r]]);
                    fdo << toRandomPeak[v[r]];
                    if(r < toRandomPeak.size()) fdo << ",";
                    below_random_threshold++;
                }
                fdo << endl;
            }
            //cout << nulldistribution.size() << endl;
        }
    } // QTL for loop

}

void analysis_cpp::functionalEnrichment(string fout){
    /*
        Count how many eQTLs are overlapping with Peak and create a 2x2 matrix that looks like this.
             | eQTL | null | 
        over |  A  |  B    |    
        no   |  C  |  D    |

        Then perform fisher test (need to code this one (: )
        */
    output_file fdo (fout);
    fdo << "A\tB\tC\tD\tOddsRatio\tpvalue" << endl;

    int qtl_overlap = 0;
    int qtl_no_overlap =0;
    int null_overlap = 0;
    int null_no_overlap=0;
    for(int i=0; i < qtl_id.size(); i++){
        if(find(intersection_snp.begin(),intersection_snp.end(), qtl_id[i]) != intersection_snp.end()){
            qtl_overlap++;
        }else{
            qtl_no_overlap++;
        }
    }

    for(int i=0; i < nulldistribution.size(); i++){
        if(find(intersection_snp.begin(), intersection_snp.end(), nulldistribution[i]) != intersection_snp.end()){
            null_overlap++;
        }else{
            null_no_overlap++;
        }
    }
    if (qtl_overlap == 0 || null_overlap == 0) {
        qtl_overlap = qtl_overlap + 1;
        null_overlap = null_overlap + 1;
    }
    double fisher = fisher_test(qtl_overlap,null_overlap, qtl_no_overlap,null_no_overlap);
    double oddsratio = odds_ratio(qtl_overlap,null_overlap, qtl_no_overlap,null_no_overlap);

    cout << "qtl overlap: " << qtl_overlap << "\nqtl no overlap: " << qtl_no_overlap << "\nnull_overlap: " << null_overlap << "\nnull_no_overlap: " << null_no_overlap << endl;
    cout << fisher << endl;
    cout << oddsratio << endl;
    std::ostringstream strs;
    strs << fisher;


    fdo << to_string(qtl_overlap) << "\t" << to_string(null_overlap) << "\t" << to_string(qtl_no_overlap) << "\t" << to_string(null_no_overlap) << "\t" << to_string(oddsratio) << "\t" << strs.str() << endl;
}

