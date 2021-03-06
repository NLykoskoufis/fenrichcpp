#include "analysis_data.h"

using namespace std;

void analysis_cpp::readIntersection(string finter){

    //PROFILE_FUNCTION();

    input_file fd (finter);
    unsigned int linecount=0;

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
    intersection_count = linecount;
}