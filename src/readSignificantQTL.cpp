
#include "fenrichcpp.h"


using namespace std;

void fenrich_cpp::readSignificantQTL(string fnom){

    cout << "Reading [ " << fnom << " ]" << endl;

    input_file fd (fnom); // input file
    unsigned int linecount =0;

    // Read significant QTLs 
    vector < string > line; 
    string buffer;
    while(getline(fd, buffer)){
        qtl_count++;
        linecount++;
        qtl_id.push_back(buffer);
}   
    qtl_count = linecount;
    cout << "Read " << to_string(qtl_count) << " significant QTLs." << endl;
    cout << "Read " << to_string(linecount) << " lines." << endl;
}