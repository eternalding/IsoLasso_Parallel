#ifndef PROCESS_SAM_HPP
#define PROCESS_SAM_HPP
#include <string>
#include "SAMrecord.hpp"
#include <utils/ArgParser.hpp>
#include "SAM_module/ReadGroup.hpp"

namespace IsoLasso::utils
{
    bool
    ReadSamFile(const IsoLasso::format::Args& arguments);

    void
    ProcessReadGroup(const IsoLasso::format::ReadGroup& RG);

}//end of namespace IsoLasso::utils



#endif