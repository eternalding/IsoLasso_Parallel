#include <iostream>
#include "SAM_module/SAMrecord.hpp"
#include "SAM_module/SAMprocess.hpp"
#include "utils/Commontype.hpp"
#include "utils/ArgParser.hpp"
#include <string>
#include "utils/Auxiliario.hpp"
#include <chrono>
#include <ctime>

int main(int argc, char* argv[]) 
{

    //Initialization
    IsoLasso::utils::ShowInitMessage();

    //Argument Parsing
    IsoLasso::format::Args arguments;
    IsoLasso::utils::ParseArg(argc,argv,arguments);
    
    //Write ReadGroup Statistics
    auto Start = IsoLasso::utils::OpenAuxFiles(argc,argv,arguments.OUTPUT_PREFIX);

    auto [RGcount,readcount] = IsoLasso::utils::ReadSamFile(arguments);

    auto End             = IsoLasso::utils::CloseAuxFiles();
    auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(End-Start);
    auto end_time        = std::chrono::system_clock::to_time_t(End);

    std::cout<<"Total ReadGroup:"<<RGcount  <<std::endl;
    std::cout<<"Total ReadCount:"<<readcount<<std::endl;
    std::cout<<"Computation finished at "   <<std::ctime(&end_time)
             <<"Elapsed time: "  << elapsed_seconds.count() << "s\n";

    return 0;
    
}