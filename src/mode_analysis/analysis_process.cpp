#include "analysis_data.h"

using namespace std;

void analysis_cpp::createNullDistribution(){

    /*
    For each eQTL:
        1. Loop over all null and get all the SNPs matching MAF, distance and that are not eQTLs (nominal != 1).
        2. From vector randomly select 10 of them.
        3. For the next eQTL, check that the eQTLs included in vector are not the same as the previous one. This step might be long maybe its better to save these ones inside a dictionary for faster search than a vector. 

    When the null distribution of SNPs are done, then we can perform the functional enrichment. 
    */

   for(int i=0; i<qtl_count; i++){
       vector < string > toRandomPeak;
       cout << "Processed " << i+1 << " eQTLs." << endl;
       for(int s=0; s < null_count; s++){
           // Check whether eQTL dist and MAF windows are matching with any TEs. 
           // If they do not match, through message
           
           if(qtl_dist_phe_var[i] > 0){
               if(qtl_dist_phe_var_from[i] <= downstream_distance[s] && qtl_dist_phe_var_to[i]>= downstream_distance[s]){
                   if(qtl_maf_from[i] <= null_maf[s] && qtl_maf_to[i] >= null_maf[s]){
                       if(nominal[s] != 1 && find(nulldistribution.begin(),nulldistribution.end(), null_id[s]) == nulldistribution.end()){
                           toRandomPeak.push_back(null_id[s]);
                       }
                   }
               }
           }else{
               if(qtl_dist_phe_var_from[i] <= upstream_distance[s] && qtl_dist_phe_var_to[i]>= upstream_distance[s]){
                    if(qtl_maf_from[i] <= null_maf[s] && qtl_maf_to[i] >= null_maf[s]){
                        if(nominal[s] != 1 && find(nulldistribution.begin(),nulldistribution.end(),null_id[s]) == nulldistribution.end()){
                            toRandomPeak.push_back(null_id[s]);
                        }
                    }
                }
            }
        }// NULL for loop
        if(toRandomPeak.size() == 0){
            cout << "No variants found" << endl;
            empty++;
            /*cout << "Increasing window_size to 3000 and MAF 0.02 and rechecking".
            for(int s=0; s < null_count; s++){
                // Check whether eQTL dist and MAF windows are matching with any TEs. 
                // If they do not match, through message
                float new_maf_from = qtl_maf[i] - (qtl_maf[i] * 0.02)
                    float new_maf_to = qtl_maf[i]
                if(qtl_dist_phe_var[i] > 0){
                    int new_qtl_dist_from = (qtl_dist_phe_var[i] - 3000) < 0 ? 0: qtl_dist_phe_var[i] - 3000;
                    int new_qtl_dist_to = qtl_dist_phe_var[i] + 3000;
                    if(new_qtl_dist_from <= downstream_distance[s] && new_qtl_dist_to >= downstream_distance[s]){
                        if(new_maf_from <= null_maf[s] && null_maf_to >= null_maf[s] && nominal[s] != 1){
                            if(ind(nulldistribution.begin(),nulldistribution.end(),null_id[i]) == nulldistribution.end()))
                        }
                    }
                }
            }
            //cout << toRandomPeak.size() << endl;*/
        }else{
            cout << toRandomPeak.size() << endl;
            std::shuffle(toRandomPeak.begin(), toRandomPeak.end(), std::default_random_engine(seed));
            if(toRandomPeak.size() >= 10){
                for(int i=0; i < 10; i++) nulldistribution.push_back(toRandomPeak[i]);
            }else{
                for(int i=0; i < toRandomPeak.size(); i++) nulldistribution.push_back(toRandomPeak[i]);
                below_random_threshold++;
            }
            cout << nulldistribution.size() << endl;
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
    for(int i=0; i < qtl_count; i++){
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
    cout << "qtl overlap: " << qtl_overlap << "\nqtl no overlap: " << qtl_no_overlap << "\nnull_overlap: " << null_overlap << "\nnull_no_overlap: " << null_no_overlap << endl;
    cout << fisher_test(qtl_overlap,null_overlap, qtl_no_overlap,null_no_overlap) << endl;
    cout << odds_ratio(qtl_overlap,null_overlap, qtl_no_overlap,null_no_overlap) << endl;

    fdo << to_string(qtl_overlap) << "\t" << to_string(null_overlap) << "\t" << to_string(qtl_no_overlap) << "\t" << to_string(null_no_overlap) << "\t" << to_string(odds_ratio(qtl_overlap,null_overlap, qtl_no_overlap,null_no_overlap)) << "\t" << to_string(fisher_test(qtl_overlap,null_overlap, qtl_no_overlap,null_no_overlap)) << endl;
}

