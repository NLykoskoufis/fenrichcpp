

#include "fenrichcpp.h"

using namespace std;


void fenrich_cpp::fenrichcpp_createTEnull(string fout){
    output_file fdo (fout); // output_file
    fdo << "#chr\tstart\tend\tid\tn_var_in_cis\tnominalQTL\tMAF\tupstream_distance\tdownstream_distance\tupstream_phen\tdownstream_phen" << endl;
    
    int cis_window=1000000;
    
    // Looping through genotypes 
    for (int t=0; t < genotype_count ; t++){
        if (t % 1000000 == 0) cout << "Processed " << to_string(t) << " variants" << endl;
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
        unsigned int nom; 
        if(nomQTL.find(genotype_id[t]) != nomQTL.end()) nom =1;
        else nom = 0;
        
        // Looping through phenotypes
        for(int g = 0; g < phenotype_count; g++){
            if (genotype_chr[t] != phenotype_chr[g]) continue;
            int tss_from = ((tss-cis_window) >0) ? (tss-cis_window) : 0;
            int tss_to = (tss+cis_window);
            int nominal;
            
            // Check whether SNP is an eQTL. If yes, jump to the next one. 
            

            if(tss_from <= phenotype_end[g] && tss_to >= phenotype_end[g]){
                id.push_back(phenotype_id[g]);
                var_in_cis++; 
                dTSS= abs(tss - phenotype_end[g]);
                if (phenotype_strand[g] == "+" && phenotype_end[g] >= tss) dTSS = - dTSS;
                if (phenotype_strand[g] == "-" && phenotype_end[g] <= tss) dTSS = - dTSS;   
                
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
        // "#chr start end id n_var_in_cis nominalQTL MAF upstream_distance downstream_distance upstream_phen downstream_phen"
        if(n_var_down == 0 && n_var_up == 0){
            fdo << genotype_chr[t] << "\t" << to_string(genotype_start[t]) << "\t" << to_string(genotype_end[t]) << "\t" << genotype_id[t] << "\t" << to_string(var_in_cis) << "\t" << to_string(nom) << "\t" << to_string(MAF) << "\t" << "NA" << "\t" << "NA" << "\t" << "NA" << "\t" << "NA" << endl;
        }else if(n_var_down != 0 && n_var_up == 0){
            fdo << genotype_chr[t] << "\t" << to_string(genotype_start[t]) << "\t" << to_string(genotype_end[t]) << "\t" << genotype_id[t] << "\t" << to_string(var_in_cis) << "\t" << to_string(nom) << "\t" << to_string(MAF) << "\t" << "NA" << "\t" << "NA" << "\t" << to_string(downstream_distance) << "\t" << downstream_phenotype << endl;
        }else if(n_var_down == 0 && n_var_up != 0){
            fdo << genotype_chr[t] << "\t" << to_string(genotype_start[t]) << "\t" << to_string(genotype_end[t]) << "\t" << genotype_id[t] << "\t" << to_string(var_in_cis) << "\t" << to_string(nom) << "\t" << to_string(MAF) << "\t" << to_string(upstream_distance) << "\t" << upstream_phenotype << "\t" << "NA" << "\t" << "NA" << "\t" << endl;
        }else{
            fdo << genotype_chr[t] << "\t" << to_string(genotype_start[t]) << "\t" << to_string(genotype_end[t]) << "\t" << genotype_id[t] << "\t" << to_string(var_in_cis) << "\t" << to_string(nom) << "\t" << to_string(MAF) << "\t" << to_string(upstream_distance) << "\t" << upstream_phenotype << "\t" << to_string(downstream_distance) << "\t" << downstream_phenotype << "\t" << endl;
        }
    }
}

int  main(int argc, char** argv){

    boost::program_options::options_description options ("Allowed options");
    options.add_options()
        ("help","produce help message")
        ("vcf", boost::program_options::value<string>(), "genotypes in VCf format.")
        ("qtl",boost::program_options::value<string>(), "Significant QTLs one SNP per line.")
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

    if (vm.count("qtl")){
        cout << "Input qtl file: " << vm["qtl"].as<string>() << endl;
    }else{
        cout << "No input file specified. You need to specify a qtl file." << endl;
    }

    if (vm.count("out")){
        cout << "Output file: " << vm["out"].as<string>() << endl;
    }else{
        cout << "You need to specify an output file." << endl;
    }

    if (vm.count("bed") && vm.count("vcf") && vm.count("out")){
        fenrich_cpp D;
        D.readPhenotypes(vm["bed"].as<string>()); // GENES 
        D.readSignificantQTL(vm["qtl"].as<string>()); // Read nominal QTLs
        D.readGenotypes(vm["vcf"].as<string>()); // TES
        D.fenrichcpp_createTEnull(vm["out"].as<string>()); // Run analysis
         
    }

    return 0;
}
