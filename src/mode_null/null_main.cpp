#include "null_data.h"



void null_main(std::vector < std::string > & argv) {
	fenrich_cpp D;

    boost::program_options::options_description opt_basic ("\x1B[32mBasics\33[0m");
	opt_basic.add_options()
		("help", "Produces option description")
		("seed", boost::program_options::value< unsigned int >()->default_value(123456), "Random number seed. Useful to replicate runs.");

    boost::program_options::options_description opt_files ("\x1B[32mI/O\33[0m");
	opt_files.add_options()
		("vcf", boost::program_options::value< std::string >(), "Genotypes in VCF/BCF format.")
		("bed", boost::program_options::value< std::string >(), "Phenotypes in BED format.")
		("qtl", boost::program_options::value< std::string >(), "QTLs in QTLtools format.")
		("out", boost::program_options::value< std::string >(), "Output file.");

    boost::program_options::options_description opt_parallel ("\x1B[32mParallelization\33[0m");
	opt_parallel.add_options()
		("region", boost::program_options::value< string >(), "Region of interest.");

    D.option_descriptions.add(opt_basic).add(opt_files).add(opt_parallel);

    //-------------------
	// 2. PARSE OPTIONS
	//-------------------
	try {
		boost::program_options::store(boost::program_options::command_line_parser(argv).options(D.option_descriptions).run(), D.options);
		boost::program_options::notify(D.options);
	} catch ( const boost::program_options::error& e ) {
		std::cerr << "Error parsing [null] command line :" << string(e.what()) << std::endl;
		exit(0);
	}
    
    std::cout << "Generating null file.";
    if(D.options.count("help")) {
        std::cout << D.option_descriptions << std::endl;
        exit(0);
    }  

    if (D.options.count("bed") && D.options.count("vcf") && D.options.count("out") && D.options.count("qtl")){
        if (D.options.count("region")) D.genomic_region = D.options["region"].as<std::string>();
        else D.genomic_region = "NA";
        D.readPhenotypes(D.options["bed"].as<std::string>()); // GENES 
        D.readSignificantQTL(D.options["qtl"].as<std::string>()); // Read nominal QTLs
        D.readGenotypes(D.options["vcf"].as<std::string>()); // TES
        D.fenrichcpp_createTEnull(D.options["out"].as<std::string>()); // Run analysis
    }

}