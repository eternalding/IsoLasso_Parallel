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

    IsoLasso::utils::ShowInitMessage();

    IsoLasso::format::Args arguments;
    IsoLasso::utils::ParseArg(argc,argv,arguments);
    
    auto Start = IsoLasso::utils::OpenAuxFiles(argc,argv,arguments.OUTPUT_PREFIX);

    auto [RGcount,readcount] = IsoLasso::utils::ReadSamFile(arguments);

    auto End = IsoLasso::utils::CloseAuxFiles();

    std::chrono::duration<double> elapsed_seconds = End-Start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(End);

    std::cout<<"Total ReadGroup:"<<RGcount<<std::endl;
    std::cout<<"Total ReadCount:"<<readcount<<std::endl;
    std::cout << "Computation finished at " << std::ctime(&end_time)
              << "Elapsed time: " << elapsed_seconds.count() << "s\n";

    return 0;
    
}