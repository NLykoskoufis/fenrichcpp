
#include "fenrichcpp.h"

using namespace std;


void fenrich_cpp::readPhenotypes(string fbed){

    input_file fd (fbed); // input file
    unsigned int linecount =0;
    unsigned int positive_strd = 0;
    unsigned int negative_strd =0;
    
    //read header 
    vector < string > str;
    string buffer;
    getline(fd, buffer);
    boost::split(str, buffer, boost::is_any_of("\t"));
    
    vector < string > line; 
    // Reading bed file line by line. 
    while(getline(fd, buffer)){
        
        boost::split(line, buffer, boost::is_any_of("\t"));
        // Check +/- negative strand 
        if(line[5] =="+"){positive_strd++;}else{negative_strd++;}

        // POPULATE PHENOTYPE VARIABLES
        phenotype_chr.push_back(line[0]);
        phenotype_start.push_back(stoi(line[1])+1);
        phenotype_end.push_back(stoi(line[2]));
        phenotype_id.push_back(line[3]);
        phenotype_strand.push_back(line[5]);
        
        phenotype_val.push_back(vector < float > (line.size()-6, 0.0));
        for (int t=6; t < line.size(); t++){
            phenotype_val.back()[t-6] = stof(line[t]);
        }
        linecount++;
    }
    phenotype_count =  phenotype_id.size();
    fd.close();    cout << to_string(phenotype_count) << endl;
    cout << positive_strd << " phenotypes on the positive strand" << endl;
    cout << negative_strd << " phenotypes on the negative strand" << endl;
}


