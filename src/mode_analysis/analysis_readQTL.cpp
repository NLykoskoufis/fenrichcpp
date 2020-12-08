#include "analysis_data.h"

using namespace std;



void analysis_cpp::readQTL(string fqtl){
    input_file fd (fqtl);
    unsigned int linecount =0;

    //read header 
    vector < string > header;
    string buffer;
    getline(fd, buffer);
    boost::split(header, buffer, boost::is_any_of(" "));
    //if(!header[0].compare("phe_id")){
    //    cout << "There is a problem with the header.Either it does not exist or it is wrong. Are you sure you are using a file generated from the QTLtools cis mode?" << endl;
    //    exit(-1);
    //}
    // CHECK TO SEE WHAT TYPE OF QTLTOOLS OUTPUT IS IT. 
    // It can be nominal, permutation, independent. The last columns will change
    
    // Read QTL information
    vector < string > line; 
    while(getline(fd, buffer)){
        if (linecount % 1000 == 0) cout << "Read " << to_string(linecount) << " lines" << endl;
        linecount++;
       
        boost::split(line, buffer, boost::is_any_of(" ")); // Split line
        qtl_id.push_back(line[7]);
        qtl_dist_phe_var.push_back(stoi(line[6]));
        
        // We need to also write the maf
        unordered_map<std::string,float>::iterator got = map_maf.find(line[7]);
        if(got == map_maf.end()){
            cout << "Problem!!! Could not find SNP. DID you use the same VCF for the creation of the null and the eQTL analysis?" << endl;
        }else{
            qtl_maf.push_back(got->second); 
            qtl_maf_from.push_back((got->second - (got->second * window_maf)));
            qtl_maf_to.push_back((got->second + (got->second * window_maf)));
        }
        
        // If d > 0 --> Downstream and if from is smaller than 0 meaning that bin is tresspassing into upstream territory, set "to" to 0. We do not want to go upstream if variant is downstream. The same for second if but for downstream.
        if(stoi(line[6]) < 0){
            if((stoi(line[6]) + window_size) > 0){
                qtl_dist_phe_var_to.push_back(0);
            }else{
                qtl_dist_phe_var_to.push_back((stoi(line[6]) + window_size));
            }
            qtl_dist_phe_var_from.push_back(stoi(line[6])-window_size);
        }else{
            if((stoi(line[6]) - window_size) < 0){
                qtl_dist_phe_var_from.push_back(0);
            }else{
                qtl_dist_phe_var_from.push_back(stoi(line[6])-window_size);
            }
            qtl_dist_phe_var_to.push_back(stoi(line[6])+window_size);
        }
    
        qtl_phe_strand.push_back(line[4]);
    }
    qtl_count = linecount;
    cout << "Read " << qtl_count << endl;
}