#ifndef PROCESS_SAM_HPP
#define PROCESS_SAM_HPP
#include <string>
#include "SAMrecord.hpp"
#include <utils/ArgParser.hpp>
#include "SAM_module/ReadGroup.hpp"

//Memory mapping
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <cstring>

namespace IsoLasso::utils
{
    std::tuple<uint32_t,uint32_t>
    ReadSamFile(const IsoLasso::format::Args& arguments);
    

    inline const char*
    MemoryMappingFile(const char* FNAME,struct stat& sb)
    { 
        auto Start        {std::chrono::system_clock::now()};
        auto start_time_t {std::chrono::system_clock::to_time_t(Start)};
        std::cout<<"Memory mapping input sam files:"<<std::endl;
        
        auto fp = open(FNAME,O_RDONLY);
        if(fp==-1)
            throw std::invalid_argument("File opening error!");
        auto FPStat = fstat(fp,&sb);
        if(FPStat==-1)
            throw std::invalid_argument("File status extraction error!");

        posix_fadvise(fp, 0, 0, POSIX_FADV_SEQUENTIAL);
        static const char* SamFile = static_cast<char*>(mmap(NULL, sb.st_size, PROT_READ, MAP_PRIVATE, fp, 0));

        if(SamFile == MAP_FAILED)
            throw std::invalid_argument("Memory mapping error!");

        auto End     {std::chrono::system_clock::now()};
        auto end_time_t {std::chrono::system_clock::to_time_t(End)};
        auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(End-Start);
        std::cout<<"Memory mapping:"<<elapsed_seconds.count()<<" sec." <<std::endl;       

        return SamFile;
    }


}//end of namespace IsoLasso::utils



#endif