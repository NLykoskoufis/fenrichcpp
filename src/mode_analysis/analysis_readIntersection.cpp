#include "analysis_data.h"

using namespace std;

void analysis_cpp::readIntersection(string finter){

    input_file fd (finter);
    unsigned int linecount;

    // Read Data 
    vector < string > line;
    string buffer;
    while(getline(fd, buffer)){
        if (linecount % 1000000 == 0) cout << "Read " << to_string(linecount) << " lines" << endl;
        linecount++;
        
        boost::split(line, buffer, boost::is_any_of("\t"));
        peaks.insert(line[9]);
        intersection_snp.push_back(line[3]);
        intersection_peaks.push_back(line[9]);
    }
    cout << "Read " << linecount << endl;
    cout << "Found " << peaks.size() << endl;

    cout << "[";
    for(auto x:peaks){
       cout << x << ","; 
    }
    cout << "] peaks." << endl;
}