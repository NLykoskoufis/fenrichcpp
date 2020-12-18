
#include "mode_null/null_data.h"
#include "mode_analysis/analysis_data.h"

void printModes(){
    cout << "Usage:" << endl;
    cout << " fenrich [mode] [options]" << endl;
    cout << " eg: fenrich null --help" << endl;
    cout << "\x1B[36;1m" << "Available modes:" << "\033[0m" <<  endl;
    cout << "\x1B[37;1m" <<   " null" << "\033[0m" << " Create null file for enrichment later on" << endl;
    cout << "\x1B[37;1m" << " enrich" << "\033[0m" << " Perform functional enrichment for QTLs" << endl;
}

int main(int argc, char ** argv) {

    cout << "\n" << "\x1B[32;1m" << "Functional enrichment for QTL variants" << "\033[0m" << endl; 
    cout << " * Authors : Nikolaos M.R. LYKOSKOUFIS" << endl;
    cout << " * Contact : nikolaos.lykoskoufis@unige.ch" << endl;
    cout << " * Version : version 1.0" << endl;

    // MODES 
    vector < string > args;
    if(argc < 2){
        printModes();
        exit(0);
    }
    for(int a= 2; a < argc ; a++) args.push_back(string(argv[a]));

    // NULL MODE 
    if (strcmp(argv[1], "null") == 0) null_main(args);

    // ENRICH MODE 
    else if (strcmp(argv[1], "enrich") == 0) analysis_main(args);

    else if (strcmp(argv[1], "--help") == 0){
        printModes();
        exit(0);
    } else {
        printModes();
        cout << "Unrecognized Fenrich mode!" << endl;
    }
}