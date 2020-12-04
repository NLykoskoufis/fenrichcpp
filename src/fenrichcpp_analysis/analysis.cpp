#include "analysis.h"

using namespace std;

int main(){
    
    string fnom;
    fnom="/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/eqtls/mapping/tumor/TUMOR.conditional_Geness_bestPerRank.txt.gz";
    string fnull="/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/eqtls/enrichment/analysis_functionalEnrichment/testing_cpp/TUMOR_variantNull_TESTING_MAP.txt.gz";
    string finter = "/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/eqtls/enrichment/analysis_functionalEnrichment/intersect/LoVo_ChIPseq_IRX3_variant_overlap.bed.gz";
    string fout = "test.out.txt";

    analysis_cpp D;
    D.readNull(fnull);
    D.readIntersection(finter);
    D.readQTL(fnom);
    D.createNullDistribution();
    D.functionalEnrichment(fout);
    /*
    cout << "Reading results loaded in memory." << endl; 
    cout << D.qtl_id.size() << endl;
    cout << D.window_maf << endl; 
    cout << D.window_size << endl;
    for (int i = 0; i < D.qtl_id.size(); i++){
        cout << D.qtl_maf[i] << " " << D.qtl_maf_from[i] << " " << D.qtl_maf_to[i] << " " << D.qtl_dist_phe_var[i] << " " << D.qtl_dist_phe_var_from[i] << " " << D.qtl_dist_phe_var_to[i] << endl;
    }*/
    //for(int i = 0; i<D.nulldistribution.size(); i++){ cout << D.nulldistribution[i] << " ";
    //} cout << endl;
    cout << D.nulldistribution.size() << endl;
    cout << D.empty << endl;
    cout << D.below_random_threshold << endl;
    
    


    return 0;
}    