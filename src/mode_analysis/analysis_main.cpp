#include "analysis_data.h"

using namespace std;

/*void analysis_main(vector < string > & argv){
    analysis_cpp D;

    boost::program_options::options_description options ("Allowed options");
    options.add_options()
        ("help","produce help message")
        ("null", boost::program_options::value<string>(), "null file created with null mode.")
        ("qtl",boost::program_options::value<string>(), "QTL output file. Same layout as QTLtools cis output.")
        ("phen", boost::program_options::value<string>(), "phenotypes to enrich. Specific format obtained from bedtools intersect [For now. Maybe later will do the overlap inside the code :)]")
        ("out", boost::program_options::value<string>(), "Output file.")
    ;

    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::parse_command_line(argv, options), vm);
    boost::program_options::notify(vm);

    if (vm.count("help")) {
        cout << options << endl;
    }

    if (vm.count("null")){
        cout << "Input null file: " << vm["vcf"].as<string>() << endl;
    }else{
        cout << "No input file specified. You need to specify a null file." << endl;
    }

    if (vm.count("phen")){
        cout << "Input phenotype file: " << vm["bed"].as<string>() << endl;
    }else{
        cout << "No input file specified. You need to specify a phenotype file." << endl;
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
        D.readNull(vm["null"].as<string>()); 
        D.readIntersection(vm["phen"].as<string>()); 
        D.readQTL(vm["qtl"].as<string>()); 
        cout << "Creating null distribution for the enrichment" << endl;
        D.createNullDistribution();
        cout << "Performing Enrichment over the null" << endl;
        D.functionalEnrichment(vm["out"].as<string>());
    }
}
*/

void analysis_main(vector < string > & argv) {
	analysis_cpp D;
    boost::program_options::options_description opt_basic ("\x1B[32mBasics\33[0m");
	opt_basic.add_options()
		("help", "Produces option description")
		("seed", boost::program_options::value< unsigned int >()->default_value(123456), "Random number seed. Useful to replicate runs.")
		("log", boost::program_options::value< string >(), "Output on screen goes to this file.")
		("silent", "Disable screen output");

    boost::program_options::options_description opt_files ("\x1B[32mI/O\33[0m");
	opt_files.add_options()
		("null", boost::program_options::value< string >(), "Null file created with fenrich null mode.")
		("phen", boost::program_options::value< string >(), "Intersection between variants and phenotypes in specific format (see documentation).")
		("qtl", boost::program_options::value< string >(), "QTLs in QTLtools format.")
		("out", boost::program_options::value< string >(), "Output file.");

	boost::program_options::options_description opt_parameters ("\x1B[32mParameters\33[0m");
    opt_parameters.add_options()
		("random_var",boost::program_options::value< unsigned int >()->default_value(10), "The number of random_variants to use")
		("window_size", boost::program_options::value< unsigned int >()->default_value(2500), "Size of the cis-window.")
        ("maf_window", boost::program_options::value< float >()->default_value(0.01), "Size of the maf_window");

    D.option_descriptions.add(opt_basic).add(opt_files).add(opt_parameters);

    //-------------------
	// 2. PARSE OPTIONS
	//-------------------
	try {
		boost::program_options::store(boost::program_options::command_line_parser(argv).options(D.option_descriptions).run(), D.options);
		boost::program_options::notify(D.options);
	} catch ( const boost::program_options::error& e ) {
		cerr << "Error parsing [enrich] command line :" << string(e.what()) << endl;
		exit(0);
	}

    std::cout << "\n" << "\x1B[32;1m" << "FUNCTIONAL ENRICHMENT ANALYSIS" << "\033[0m" << endl;
    if(D.options.count("help")) {
        cout << D.option_descriptions << endl;
        exit(0);
    }

    if (D.options.count("phen") && D.options.count("null") && D.options.count("out") && D.options.count("qtl")){
        D.random_variants = D.options["random_var"].as < unsigned int >();
        D.window_size = D.options["window_size"].as < unsigned int >();
        D.window_maf = D.options["maf_window"].as < float >();

        cout << "performing functional enrichment" << endl;
        cout << " * Random variants: " << D.random_variants << endl;
        cout << " * Window size    :" << D.window_size << endl;
        cout << " * MAF window     :" << D.window_maf << endl;

        cout << "Reading [" << D.options["null"].as<string>() << "]" << std::endl;
        D.readNull(D.options["null"].as<string>()); 

        //cout << "Reading [" << D.options["phen"].as<string>() << "]"; 
        //D.readIntersection(D.options["phen"].as<string>()); 
        
        cout << "Reading [" << D.options["qtl"].as<string>() << "]" << std::endl;
        D.readQTL(D.options["qtl"].as<string>()); 
        
        cout << "Creating null distribution for the enrichment" << endl;
        D.createNullDistribution();
        //cout << "Performing Enrichment over the null" << endl;
        //D.functionalEnrichment(D.options["out"].as<string>());
        cout << "Reading [" << D.options["phen"].as<string>() << "]" << std::endl;
        D.readPhenotypes(D.options["phen"].as<string>());
        
        std::cout << D.phen_chr[0] << std::endl;
        std::cout << D.qtl_chr[0] << std::endl;
        std::cout << D.qtl_id.size() << std::endl;
        
        std::cout << "Performing qtl intersection" << std::endl;
        D.performIntersect();
    }
}