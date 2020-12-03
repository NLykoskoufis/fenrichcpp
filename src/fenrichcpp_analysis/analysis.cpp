#include "analysis.h"

using namespace std;

int main(){
    
    string fnom;
    fnom="/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/eqtls/permute/tumor/TUMOR.permute_All.significant.txt";
    string fnull="/Users/nikolaoslykoskoufis/Documents/PROJECTS/fenrichcpp/testing_code.txt.gz";
    analysis_cpp D;
    D.readNull(fnull);
    //D.readQTL(fnom);
    return 0;
}    