#ifndef AUXILIARIO_HPP
#define AUXILIARIO_HPP
#include <SAM_module/ReadGroup.hpp>
#include <utils/Commontype.hpp>
#include <string>
#include <filesystem>
#include <chrono>
#include <iomanip>

namespace fs = std::filesystem;

namespace IsoLasso::utils
{
    inline void
    ShowInitMessage()
    {
        std::cout<<"* IsoLasso Parallel, Version:"<<version<<std::endl;
        std::cout<<"* Author: Yu-Cheng Lo"<<std::endl;
        return;
    }

    inline std::chrono::_V2::system_clock::time_point
    OpenAuxFiles(int argc,char* argv[],const std::string& OUTPUT_PREFIX)
    {
        auto Start     {std::chrono::system_clock::now()};
        auto in_time_t {std::chrono::system_clock::to_time_t(Start)};

        GTF_FS.open(OUTPUT_PREFIX+"_isoforms.gtf");
        if(!GTF_FS.is_open())
            std::__throw_invalid_argument("Cannot create .gtf file!");

        RG_STATS_FS.open(OUTPUT_PREFIX+"_RGStats.txt");
        if(!RG_STATS_FS.is_open())
            std::__throw_invalid_argument("Cannot write into ReadGroup Stats File!");
        RG_STATS_FS<<"* IsoLasso Parallel, Version:"<<version<<std::endl;
        RG_STATS_FS<<"* Start time:"<<std::put_time(std::localtime(&in_time_t), "%Y-%m-%d %X")<<std::endl;
        RG_STATS_FS<<"* Command:";
        for(auto arg_idx=0;arg_idx<argc;arg_idx++)  
            RG_STATS_FS<<argv[arg_idx]<<" ";
        RG_STATS_FS<<std::endl;
        return Start;
    }

    inline std::chrono::_V2::system_clock::time_point
    CloseAuxFiles()
    {
        auto End {std::chrono::system_clock::now()};
        auto in_time_t = std::chrono::system_clock::to_time_t(End);

        RG_STATS_FS.close();
        return End;
    }

    template <class T>
    inline void 
    print2Dvector(const std::vector<std::vector<T>>& target)
    {
        for(const auto& row:target)
        {
            for(const auto col:row)
                std::cout<<col<<" ";
            std::cout<<std::endl;
        }
        return;
    }

    template <class T>
    inline void 
    print1Dvector(const std::vector<T>& target)
    {
        for(const auto elem:target)        
            std::cout<<elem<<" ";
        std::cout<<std::endl;

        return;
    }

}
#endif