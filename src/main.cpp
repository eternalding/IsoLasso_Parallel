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
    
    IsoLasso::utils::ReadSamFile(arguments);
    //ReadSamFile(arguments.SAMFILE);

    return 0;
    
}