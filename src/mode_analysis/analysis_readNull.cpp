#include "analysis_data.h"

using namespace std;



void analysis_cpp::readNull(string fnull){
    
    PROFILE_FUNCTION();
    
    input_file fd (fnull);
    unsigned int linecount =0;

    //read header 
    vector < string > header;
    string buffer;
    getline(fd, buffer);
    boost::split(header, buffer, boost::is_any_of("\t"));
    /*if(!header[0].compare("#chr")){
        cout << "There is a problem with the header.Either it does not exist or it is wrong. " << endl;
        exit(-1);
    }*/
    
    // Read Null 
    vector < string > line; 
    while(getline(fd, buffer)){
        if (linecount % 1000000 == 0) cout << "Read " << to_string(linecount) << " lines" << endl;
        linecount++;

        boost::split(line, buffer, boost::is_any_of("\t"));
        null_id.push_back(line[3]);
        
        if(has_chr(line[0])) null_chr.push_back(line[0].substr(3,line[0].size()));
        else null_chr.push_back(line[0]);

        null_start.push_back(std::stoi(line[1]));
        null_end.push_back(std::stoi(line[2]));

        nominal.push_back(stoi(line[5]));
        null_maf.push_back(stof(line[6]));
        if(line[7].find("NA") != string::npos) upstream_distance.push_back(stoi("1"));
        else upstream_distance.push_back(stoi(line[7]));
        if(line[9].find("NA") != string::npos) downstream_distance.push_back(stoi("-1"));
        else downstream_distance.push_back(stoi(line[9]));
        //cout << line[3] << endl;
        map_maf.insert(make_pair(line[3],stof(line[6]))); // Create dictionary with SNP::MAF

    }
    null_count = linecount;
    cout << "Read " << to_string(linecount) << " null variants." << endl;
}
