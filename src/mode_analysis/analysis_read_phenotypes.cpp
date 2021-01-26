#include "analysis_data.h"

void analysis_cpp::readPhenotypes(std::string fbed)
{
    input_file fd (fbed);
    std::string buffer;
    std::vector < std::string > line;
    unsigned int linecount = 0;

    while(getline(fd, buffer)){
        if (linecount % 1000000 == 0) cout << "Read " << to_string(linecount) << " lines" << endl;

        boost::split(line,buffer, boost::is_any_of("\t"));
        linecount++;

        if(line.size() < 3)
        {
            std::cout << "Less than 3 columns. Please use BED format";
            exit(-1);
        }
        if(has_chr(line[0]))
        {
            phen_chr.push_back(line[0].substr(3,line[0].size()));
        }else{
            phen_chr.push_back(line[0]);
        }
        phen_start.push_back(std::stoi(line[1])+1); // Set to 1 based.
        phen_end.push_back(stoi(line[2]));

    }
    phen_count = linecount;
    std::cout << "Read " << phen_count << std::endl;
}