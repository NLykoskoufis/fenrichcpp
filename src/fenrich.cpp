
#include "mode_null/null_data.h"
#include "mode_analysis/analysis_data.h"

void printModes(){
    cout << "Usage:" << endl;
    cout << " fenrich [mode] [options]" << endl;
    cout << " eg: fenrich null --help" << endl;
    cout << "Available modes:" << endl;
    cout << " null  Create null file for enrichment later on" << endl;
    cout << " enrich    Perform functional enrichment for QTLs" << endl;
}

int main(int argc, char ** argv) {

    cout << "Fenrich" << endl; 
    cout << "Authors : Nikolaos M.R. LYKOSKOUFIS" << endl;
    cout << "Contact : nikolaos.lykoskoufis@unige.ch" << endl;
    cout << "Version : version 0.1" << endl;

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