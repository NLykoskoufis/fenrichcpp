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
  
            }

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
        if(find(qtl_id.begin(),qtl_id.end(), snp_id[i]) == qtl_id.end()){
        
            qtl_id.push_back(snp_id[i]);
            qtl_dist_phe_var.push_back(dist);

            // GET MAF 
            unordered_map<std::string,float>::iterator got = map_maf.find(snp_id[i]);
            if(got == map_maf.end()){
                cout << "Problem!!! Could not find SNP. DID you use the same VCF for the creation of the null and the eQTL analysis?" << endl;
            }else{
                qtl_maf.push_back(got->second); 
                }
            }
        }
    }
    cout << snp_id.size() << " " << qtl_id.size() << endl;
    cout << multiple << endl;

}

