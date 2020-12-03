#include "analysis.h"

using namespace std;



void analysis_cpp::readNull(string fnull){
    input_file fd (fnull);
    unsigned int linecount =0;

    //read header 
    vector < string > header;
    string buffer;
    getline(fd, buffer);
    boost::split(header, buffer, boost::is_any_of(" "));
    if(!header[0].compare("#chr")){
        cout << "There is a problem with the header.Either it does not exist or it is wrong. " << endl;
        exit(-1);
    }
    
    
    // Read Null 
    vector < string > line; 
    while(getline(fd, buffer)){
        if (linecount % 1000 == 0) cout << "Read " << to_string(linecount) << " lines" << endl;
        linecount++;
        null_count++;

        boost::split(line, buffer, boost::is_any_of("\t"));
        null_id.push_back(line[0]);
        null_maf.push_back(stof(line[5]));
        upstream_distance.push_back(stoi(line[6]));
        downstream_distance.push_back(stoi(line[8]));
    }   

}
