#include "analysis_data.h"

void analysis_cpp::readPhenotypes(std::string fbed)
{
    input_file fd (fbed);
    std::string buffer;
    std::vector < std::string > line;
    int linecount = 0;
    std::string chrom;
    chrom = "1";
    int index_start = 0;
    int index_end = 0;

    while(getline(fd, buffer)){
        if (linecount % 1000000 == 0) cout << "Read " << to_string(linecount) << " lines" << endl;
        linecount++;

        boost::split(line,buffer, boost::is_any_of("\t"));
        

        if(line.size() < 3)
        {
            std::cout << "Less than 3 columns. Please use BED format";
            exit(-1);
        }
        std::string chr;
        if(has_chr(line[0]))
        {
            phen_chr.push_back(line[0].substr(3,line[0].size()));
            chr = line[0].substr(3,line[0].size());
        }else{
            phen_chr.push_back(line[0]);
            chr = line[0];
        }
        phen_start.push_back(std::stoi(line[1])+1); // Set to 1 based.
        phen_end.push_back(stoi(line[2]));

        // Creating index
        if (chrom.compare(chr) == 0)
        {
            index_end = linecount;
        }else{
            if(phen_index.count(chrom))
            {
                std::cout << "Chromosome already present in index map. The file might not be properly sorted. Please sort the file using either unix sort or bedtools sort.";
                exit(-1);
            }
            positions pos = {index_start, index_end};
            phen_index.insert(std::make_pair(chrom,pos));

            chrom = chr;
            index_start = index_end + 1;
        }
    }
    positions pos = {index_start, index_end};
    phen_index.insert(std::make_pair(chrom,pos));

    phen_count = linecount;
    std::cout << "Read " << phen_count << std::endl;
    /*for(auto x:phen_index)
    {
        std::cout << x.first << " -> " << x.second.start << " " << x.second.end << std::endl;
    }*/
}