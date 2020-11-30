

#include "fenrichcpp.h"

using namespace std;

void fenrich_cpp::readGenotypes(string fvcf){
    // Opening files
    bcf_srs_t * sr = bcf_sr_init();
    if(!(bcf_sr_add_reader (sr, fvcf.c_str()))) {
		switch (sr->errnum) {
		case not_bgzf: cout << "File not compressed with bgzip!" << endl; break;
		case idx_load_failed: cout << "Impossible to load index file!"<< endl; break;
		case file_type_error: cout << "File format not detected by htslib!" << endl; break;
		default : cout << "Unknown error!" << endl;
		}
	}
    int n_samples = bcf_hdr_nsamples(sr->readers[0].header);

    unsigned int linecount=0;
    // Read genotype data 
    int ngt, ngt_arr = 0, nds, nds_arr = 0, * gt_arr = NULL, nsl, nsl_arr = 0, * sl_arr = NULL;
    float * ds_arr = NULL;
    bcf1_t * line;
    while(bcf_sr_next_line (sr)){
        linecount ++;
        genotype_count++;
        if (linecount % 100000 == 0) cout << "Read " << to_string(linecount) << " lines" << endl;
        line = bcf_sr_get_line(sr,0);
        if(line->n_allele == 2){
            ngt = bcf_get_genotypes(sr->readers[0].header, line, &gt_arr, &ngt_arr);
            if(ngt == 2*n_samples){
                string sid = string(line->d.id);
                string chr = string(bcf_hdr_id2name(sr->readers[0].header, line->rid));
                int pos = line->pos +1;

                genotype_id.push_back(sid); // Read variant ID
                genotype_chr.push_back(chr); // Read variant chr
                string genotype_ref = string(line->d.allele[0]); // Read reference allele.
                genotype_start.push_back(pos);
                nsl = bcf_get_info_int32(sr->readers[0].header, line, "END", &sl_arr, &nsl_arr);
                if (nsl >= 0 && nsl_arr == 1) genotype_end.push_back(sl_arr[0]);
                else genotype_end.push_back(genotype_start.back() + genotype_ref.size() - 1);
                genotype_val.push_back(vector < float > (n_samples, 0.0));

                for(int i = 0; i < n_samples ; i ++) {
                    if (gt_arr[2*i+0] == bcf_gt_missing || gt_arr[2*i+1] == bcf_gt_missing) bcf_float_set_missing(genotype_val.back()[i]);
                    else genotype_val.back()[i] = bcf_gt_allele(gt_arr[2*i+0]) + bcf_gt_allele(gt_arr[2*i+1]);
				}
		    } 
        }
    }
}


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
    phenotype_count = 0;
    // Reading bed file line by line. 
    while(getline(fd, buffer)){
        
        boost::split(line, buffer, boost::is_any_of("\t"));
        // Check +/- negative strand 
        if(line[5] =="+"){positive_strd++;}else{negative_strd++;}

        // POPULATE PHENOTYPE VARIABLES
        phenotype_chr.push_back(line[0]);
        phenotype_start.push_back(stoi(line[1]));
        phenotype_end.push_back(stoi(line[2]));
        phenotype_id.push_back(line[3]);
        phenotype_strand.push_back(line[5]);
        
        phenotype_val.push_back(vector < float > (line.size()-6, 0.0));
        for (int t=6; t < line.size(); t++){
            phenotype_val.back()[t-6] = stof(line[t]);
        }
        phenotype_count++;
        linecount++;
    }
    fd.close();    cout << linecount << endl;
    cout << positive_strd << " phenotypes on the positive strand" << endl;
    cout << negative_strd << " phenotypes on the negative strand" << endl;
}



void fenrich_cpp::fenrichcpp_createTEnull(string fout){
    output_file fdo (fout); // output_file
    fdo << "#chr\tstart\tend\tid\tn_var_in_cis\tMAF\tupstream_distance\tdownstream_distance\tupstream_phen\tdownstream_phen" << endl;
    
    int cis_window=1000000;
    
    // Looping through genotypes 
    for (int t=0; t < genotype_count ; t++){
        if (t % 100000 == 0) cout << "Processed " << to_string(t) << " variants" << endl;
        unsigned int var_in_cis=0;
        unsigned int n_var_down = 0;
        unsigned int n_var_up = 0;
        // get first Variant id, and TSS. 
        int tss = genotype_end[t];
        vector <string> id;
        int downstream_distance=10000000;
        int upstream_distance=-10000000;
        string downstream_phenotype;
        string upstream_phenotype;
        int dTSS;
        float MAF = getMAF(genotype_val[t]);

        // Looping through phenotypes
        for(int g = 0; g < phenotype_count; g++){
            if (genotype_chr[t] != phenotype_chr[g]) continue;
            int tss_from = ((tss-cis_window) >0) ? (tss-cis_window) : 0;
            int tss_to = (tss+cis_window);
            
            if(tss_from <= phenotype_start[g] && tss_to >= phenotype_start[g]){
                id.push_back(phenotype_id[g]);
                var_in_cis++; 
                dTSS= abs(tss - phenotype_start[g]);
                //cout << phenotype_strand[g] << " " << phenotype_start[g] << endl;
                if (phenotype_strand[g] == "+" && phenotype_start[g] >= tss) dTSS = - dTSS;
                if (phenotype_strand[g] == "-" && phenotype_start[g] <= tss) dTSS = - dTSS;   
                
                if (dTSS > 0 ){
                    if(dTSS < downstream_distance){
                        downstream_distance = dTSS;
                        downstream_phenotype = phenotype_id[g];
                        n_var_down ++;
                    }
                }else{
                    if(dTSS > upstream_distance){
                        upstream_distance = dTSS;
                        upstream_phenotype = phenotype_id[g];
                        n_var_up++;
                    }
                }
            }
        }

        if(n_var_down == 0 && n_var_up == 0){
            fdo << genotype_chr[t] << "\t" << to_string(genotype_start[t]) << "\t" << to_string(genotype_end[t]) << "\t" << genotype_id[t] << "\t" << to_string(var_in_cis) << "\t" << to_string(MAF) << "\t" << "NA" << "\t" << "NA" << "\t" << "NA" << "\t" << "NA" << endl;
        }else if(n_var_down != 0 && n_var_up == 0){
            fdo << genotype_chr[t] << "\t" << to_string(genotype_start[t]) << "\t" << to_string(genotype_end[t]) << "\t" << genotype_id[t] << "\t" << to_string(var_in_cis) << "\t" << to_string(MAF) << "\t" << to_string(upstream_distance) << "\t" << upstream_phenotype << "\t" << "NA" << "\t" << "NA" << endl;
        }else if(n_var_down == 0 && n_var_up != 0){
            fdo << genotype_chr[t] << "\t" << to_string(genotype_start[t]) << "\t" << to_string(genotype_end[t]) << "\t" << genotype_id[t] << "\t" << to_string(var_in_cis) << "\t" << to_string(MAF) << "\t" << "NA" << "\t" << "NA" << "\t" << to_string(downstream_distance) << "\t" << downstream_phenotype << "\t" << endl;
        }else{
            fdo << genotype_chr[t] << "\t" << to_string(genotype_start[t]) << "\t" << to_string(genotype_end[t]) << "\t" << genotype_id[t] << "\t" << to_string(var_in_cis) << "\t" << to_string(MAF) << "\t" << to_string(upstream_distance) << "\t" << upstream_phenotype << "\t" << to_string(downstream_distance) << "\t" << downstream_phenotype << "\t" << endl;
        }
    }
}

int  main(int argc, char** argv){

    boost::program_options::options_description options ("Allowed options");
    options.add_options()
        ("help","produce help message")
        ("vcf", boost::program_options::value<string>(), "genotypes in VCf format.")
        ("bed", boost::program_options::value<string>(), "phenotypes in BED format.")
        ("out", boost::program_options::value<string>(), "Output file.")
    ;

    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, options), vm);
    boost::program_options::notify(vm);

    if (vm.count("help")) {
        cout << options << endl;
    }

    if (vm.count("vcf")){
        cout << "Input vcf file: " << vm["vcf"].as<string>() << endl;
    }else{
        cout << "No input file specified. You need to specify a vcf file." << endl;
    }

    if (vm.count("bed")){
        cout << "Input bed file: " << vm["bed"].as<string>() << endl;
    }else{
        cout << "No input file specified. You need to specify a bed file." << endl;
    }

    if (vm.count("out")){
        cout << "Output file: " << vm["out"].as<string>() << endl;
    }else{
        cout << "You need to specify an output file." << endl;
    }

    if (vm.count("bed") && vm.count("vcf") && vm.count("out")){
        fenrich_cpp D;
        D.readPhenotypes(vm["bed"].as<string>()); // GENES 
        D.readGenotypes(vm["vcf"].as<string>()); // TES
        D.fenrichcpp_createTEnull(vm["out"].as<string>());
    }

    return 0;
}
