#include <SAM_module/SAMrecord.hpp>
#include <SAM_module/SAMprocess.hpp>
#include <SAM_module/ReadGroup.hpp>
#include <utils/Commontype.hpp>
#include <utils/ArgParser.hpp>
#include <utils/Auxiliario.hpp>
#include <iostream>
#include <fstream>
#include <numeric>
#include <utils/ThreadPool.hpp>
#include <thread>
#include <string_view>

//Memory mapping
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <cstring>
#include <streambuf>
#include <istream>
#include <sstream>
#include <string_view>
#include <iostream>

namespace IsoLasso::utils
{
    std::vector<std::string_view>
    splitSV(std::string_view strv, std::string_view delims = " ")
    {
        std::vector<std::string_view> output;
        size_t first = 0;

        while (first < strv.size())
        {
            const auto second = strv.find_first_of(delims, first);

            if (first != second)
                output.emplace_back(strv.substr(first, second-first));

            if (second == std::string_view::npos)
                break;

            first = second + 1;
        }

        return output;
    }


    std::tuple<uint32_t,uint32_t>
    ReadSamFile(const IsoLasso::format::Args& arguments)
    {
        /*        
        struct stat sb;

        auto MMFStart = std::chrono::system_clock::now();
        auto SamFile  = IsoLasso::utils::MemoryMappingFile(arguments.SAMFILE.c_str(),sb);
        auto TotalSize = SamFile+sb.st_size;
        std::cout<<sb.st_size<<std::endl;

        std::stringstream iss;

        //std::cout<<SamFile<<std::endl;
        //std::cout<<*TotalSize<<std::endl;

        std::cout<<"TotalSize:"<<TotalSize<<std::endl;

        auto reads = splitSV(std::string_view(SamFile),"\n");

        std::cout<<"LineCount:"<<reads.size()<<std::endl;

        std::cout<<reads[2000]<<std::endl;
        std::cout<<reads[8000]<<std::endl;

        auto Start = std::chrono::system_clock::now();
        */
        std::ifstream fin(arguments.SAMFILE,std::ios::in);
        //std::ios::sync_with_stdio(false);
        if (!fin.is_open())
            throw std::invalid_argument("Failed to open target SAM file!");

        IsoLasso::format::Header_record HRecord;
        IsoLasso::format::Sam_record Record;
        IsoLasso::format::ReadGroup RG;
        uint32_t Total_RG_index  {0};
        uint32_t Valid_RG_index  {0};
        uint32_t ReadCount       {0};
        std::vector<uint32_t> rg_size;

        //ThreadPool
        const auto processor_count = std::thread::hardware_concurrency();
        thread_pool pool(processor_count);

        //Header is now unavailable
        while(fin>>HRecord);
        
        //Main workflow
        while(fin>>Record)
        {
            //SAM file should be sorted
            if(Record.Pos<RG.CurrentRange.first)
            {
                std::cerr<<"Exception: "<<Record.Pos<<"/"<<RG.CurrentRange.first<<std::endl;
                std::__throw_invalid_argument("Using unsorted SAM file! Aborting...");
            }

            // Ignore invalid SAM records
            if(Record.ValidBit)
            {
                //Ignore empty RG 
                if(RG.size()==0)
                {}
                //Begin new readGroup
                else if(Record.RName != RG.ChrName || Record.Pos > RG.CurrentRange.second + MIN_GAP_SPAN)
                {                
                    RG.RG_index = Total_RG_index;
                    if(RG.size()>MIN_RG_SIZE)
                    {
#ifdef DEBUG
                        std::cout<<"Processing ReadGroup:"<<RG.RG_index
                                 <<", Range:"<<RG.ChrName
                                 <<":["<<RG.CurrentRange.first
                                 <<","<<RG.CurrentRange.second
                                 <<"]\n";
#endif
                        Valid_RG_index++;
                        ReadCount += RG.validSize();
                        pool.submit(IsoLasso::utils::ProcessReadGroup,RG);
                    }
                    else if (RG.size()>0)// RG is not large enough
                    {
#ifdef DEBUG
                        std::cout<<"Read Group "<<RG.RG_index
                                 <<" with size "<<RG.size()
                                 <<" is not large enough!"<<std::endl;
#endif
                    }
                    Total_RG_index++;               
                    RG.reset();                    
                }

                //Initialize ReadGroup
                if(RG.ChrName=="")
                    RG.ChrName = Record.RName;
                
                //Add current record to current ReadGroup
                RG.AddRecord(Record);

                auto [current_start,current_end] = Record.GetRange();

                if(RG.CurrentRange.first==0)
                    RG.CurrentRange.first = current_start;

                RG.CurrentRange.second = (current_end>RG.CurrentRange.second)
                                         ?current_end
                                         :RG.CurrentRange.second;
#ifdef DEBUG
                if(RG.validSize()<10)
                    std::cout<<"["<<RG.CurrentRange.first<<","<<RG.CurrentRange.second<<"]"<<"\n";
#endif
            }
        }//End of reading SamFILE

        //Last ReadGroup
        RG.RG_index = Total_RG_index;
        Valid_RG_index++;

        if(RG.ReadStart.size()>MIN_RG_SIZE)
        {
            ReadCount += RG.validSize();
#ifdef DEBUG
            std::cout<<"Processing Last ReadGroup "<<RG.RG_index<<", Range:"<<RG.ChrName<<" ["<<RG.CurrentRange.first<<","<<RG.CurrentRange.second<<"]"<<std::endl;
            std::cout<<"Number of valid reads:"<<RG.validSize()<<"\n";
#endif
            ProcessReadGroup(RG);   
        }

        RG_STATS_FS.close();
        GTF_FS.close();

        return {Valid_RG_index,ReadCount};
    }//end of ReadSamFile
};//end of namespace IsoLaso::utils