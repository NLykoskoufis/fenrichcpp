#include "analysis.h"

using namespace std;

int main(){
    
    string fnom;
    fnom="/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/eqtls/permute/tumor/TUMOR.permute_All.significant.txt";
    
    analysis_cpp D;
    D.readQTL(fnom);

    }    
}