#ifndef ARGPARSER_HPP
#define ARGPARSER_HPP

#include <iostream>
#include <array>
#include "Commontype.hpp"

namespace IsoLasso::format
{
    struct Args
    {
        //Required fields
        std::string SAMFILE              ;

        //Optional fields
        std::string OUTPUT_PREFIX        ; // Specify the prefix of all generated files. 
        
    };
    extern IsoLasso::format::Args arguments;

}//end of namespace IsoLasso::format

namespace IsoLasso::utils
{
    void 
    printHelpMsg();
    
    bool 
    ParseArg(const int& argc ,char* argv[],IsoLasso::format::Args& Args);
}//end of namespace IsoLasso::utils



#endif