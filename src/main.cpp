#include <iostream>
#include "SAM_module/SAMrecord.hpp"
#include "SAM_module/SAMprocess.hpp"
#include "utils/Commontype.hpp"
#include "utils/ArgParser.hpp"
#include <string>

int main(int argc, char* argv[])
{

    IsoLasso::format::Args arguments;

    IsoLasso::utils::ParseArg(argc,argv,arguments);
    
    auto [RGcount,readcount] = IsoLasso::utils::ReadSamFile(arguments);
    
    std::cout<<"Total ReadGroup:"<<RGcount<<std::endl;
    std::cout<<"Total ReadCount:"<<readcount<<std::endl;

    return 0;
    
}