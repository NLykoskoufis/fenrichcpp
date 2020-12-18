#include "analysis_data.h"
#include <cstdlib>
using namespace std;



void analysis_cpp::readQTL(string fqtl){
    input_file fd (fqtl);
    unsigned int linecount =0;
    qtl_count = 0;

    vector < string > snp_id;
    vector < int > dist_phe_var;

    //read header 
    vector < string > header;
    string buffer;
    getline(fd, buffer);
    boost::split(header, buffer, boost::is_any_of(" "));
    
    // Read QTL information
    vector < string > line;
    while(getline(fd, buffer)){
        if (linecount % 1000 == 0) cout << "Read " << to_string(linecount) << " lines" << endl;
        linecount++;
        boost::split(line, buffer, boost::is_any_of(" ")); // Split line
        qtl_count++;
        
        snp_id.push_back(line[7]);
        dist_phe_var.push_back(stoi(line[6]));
        
    } // Finished reading the file. 

    unordered_map<string,unsigned int> qtl_check;


    cout << "Read " << linecount << endl;
    cout << "Processing QTLs. When multiple SNPs, keeping only the one closest to the gene.";
    int multiple=0;
    for(int i =0; i < linecount; i++){
        int dist = 1000000000;
        if(count(snp_id.begin(),snp_id.end(), snp_id[i])==1){
            //cout << snp_id[i] << " " << count(snp_id.begin(),snp_id.end(), snp_id[i]) << endl;
            qtl_id.push_back(snp_id[i]);
            qtl_dist_phe_var.push_back(dist_phe_var[i]);
            // We need to also write the maf
            unordered_map<std::string,float>::iterator got = map_maf.find(snp_id[i]);
            if(got == map_maf.end()){
                cout << "Problem!!! Could not find SNP. DID you use the same VCF for the creation of the null and the eQTL analysis?" << endl;
            }else{
                qtl_maf.push_back(got->second); 
                //qtl_maf_from.push_back((got->second - (got->second * window_maf)));
                //qtl_maf_to.push_back((got->second + (got->second * window_maf)));
            }

            // GET UPSTREAM AND DOWNSTREAM 
            // If d > 0 --> Downstream and if from is smaller than 0 meaning that bin is tresspassing into upstream territory, set "to" to 0. We do not want to go upstream if variant is downstream. The same for second if but for downstream.
            /*if(dist_phe_var[i] < 0){
                if((dist_phe_var[i] + window_size) > 0){
                    qtl_dist_phe_var_to.push_back(0);
                }else{
                    qtl_dist_phe_var_to.push_back((dist_phe_var[i] + window_size));
                }
                qtl_dist_phe_var_from.push_back(dist_phe_var[i]-window_size);
            }else{
                if((dist_phe_var[i]- window_size) < 0){
                    qtl_dist_phe_var_from.push_back(0);
                }else{
                    qtl_dist_phe_var_from.push_back(dist_phe_var[i]-window_size);
                }
                qtl_dist_phe_var_to.push_back(dist_phe_var[i]+window_size);
            }*/

        }else{
            multiple++;
            cout << snp_id[i] << " " << count(snp_id.begin(),snp_id.end(), snp_id[i]) << endl;
            string snp = snp_id[i];
            for(int it=0; it < linecount; it++){
                if(snp.compare(snp_id[it]) == 0){
                    if(dist < abs(dist_phe_var[it])){
                        continue;

                    }else{
                        dist = dist_phe_var[it];
                    }
                }
            }
        cout << snp_id[i] << " " << dist << endl;
        if(find(qtl_id.begin(),qtl_id.end(), snp_id[i]) == qtl_id.end()){ // This is going to take more and more time to run because the vector will incrementally increase. I will try with a unordered map and check whether the key is in or not. Might be faster. 
        //unordered_map<std::string,unsigned int>::iterator ch = qtl_check.find(snp_id[i]);
        //if(ch == qtl_check.end()){
        
            qtl_id.push_back(snp_id[i]);
            qtl_dist_phe_var.push_back(dist);

            // GET MAF 
        unordered_map<std::string,float>::iterator got = map_maf.find(snp_id[i]);
            if(got == map_maf.end()){
                cout << "Problem!!! Could not find SNP. DID you use the same VCF for the creation of the null and the eQTL analysis?" << endl;
            }else{
                qtl_maf.push_back(got->second); 
                //qtl_maf_from.push_back((got->second - (got->second * window_maf)));
                //qtl_maf_to.push_back((got->second + (got->second * window_maf)));
            }
            // GET UPSTREAM AND DOWNSTREAM 
            // If d > 0 --> Downstream and if from is smaller than 0 meaning that bin is tresspassing into upstream territory, set "to" to 0. We do not want to go upstream if variant is downstream. The same for second if but for downstream.
            /*if(dist < 0){
                if((dist + window_size) > 0){
                    qtl_dist_phe_var_to.push_back(0);
                }else{
                    qtl_dist_phe_var_to.push_back((dist + window_size));
                }
                qtl_dist_phe_var_from.push_back(dist-window_size);
            }else{
                if((dist - window_size) < 0){
                    qtl_dist_phe_var_from.push_back(0);
                }else{
                    qtl_dist_phe_var_from.push_back(dist-window_size);
                }
                qtl_dist_phe_var_to.push_back(dist+window_size);
            }*/

            }
        }
        
    }
    cout << snp_id.size() << " " << qtl_id.size() << endl;
    cout << multiple << endl;

    /*cout << "For debugging purposes, writing to file." << endl;
    output_file fdo ("windows_created.txt");
    fdo << "ID\tdist_phe_var\tMAF\tfrom\tto\tmaf_from\tmaf_to" << endl;
    for (int i =0;i < qtl_id.size(); i++)
    {
        fdo << qtl_id[i] << " " << to_string(qtl_dist_phe_var[i]) << " " << to_string(qtl_maf[i]) << " " << to_string(qtl_dist_phe_var_from[i]) << " " << to_string(qtl_dist_phe_var_to[i]) << " " << to_string(qtl_maf_from[i]) << " " << to_string(qtl_maf_to[i]) << endl; 
    }*/

}

