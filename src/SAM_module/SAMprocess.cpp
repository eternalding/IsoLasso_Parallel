//IsoLasso
#include <SAM_module/SAMrecord.hpp>
#include <SAM_module/SAMprocess.hpp>
#include <SAM_module/ReadGroup.hpp>
#include <utils/Commontype.hpp>
#include <utils/ArgParser.hpp>
#include <utils/Auxiliario.hpp>
//STL
#include <iostream>
#include <fstream>
#include <numeric>
#include <utils/ThreadPool.hpp>
#include <thread>
#include <string_view>
#include <algorithm>
#include <ctime>

//Thread number for the entire process
#ifdef DEBUG
    const auto processor_count = 1;
#else
    const auto processor_count = std::thread::hardware_concurrency();
#endif

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
        std::ifstream fin(arguments.SAMFILE,std::ios::in);
        if (!fin.is_open())
            throw std::invalid_argument("Failed to open target SAM file!");

        //Calculate total read count for RPKM calculation
        std::cout<<"Now calculating number of total reads:"<<std::endl;
        IsoLasso::utils::TOTAL_READ_CNT = IsoLasso::utils::TotalReadCnt(arguments.SAMFILE);
        std::cout<<"Total reads:"<<IsoLasso::utils::TOTAL_READ_CNT<<std::endl;

        IsoLasso::format::Header_record HRecord;
        IsoLasso::format::Sam_record Record;
        IsoLasso::format::ReadGroup RG;
        uint32_t Total_RG_index  {0};
        uint32_t Valid_RG_index  {0};
        uint32_t ReadCount       {0};
        std::vector<uint32_t> rg_size;

        //ThreadPool
        thread_pool pool(processor_count);

        //Header is now unavailable
        while(fin>>HRecord);
        
        auto linecount = 0;
        //Main workflow
        while(fin>>Record)
        {
            //SAM file should be sorted
            if(Record.RName == RG.ChrName && Record.Pos<RG.CurrentRange.first)
            {
                std::cout<<linecount<<std::endl;
                std::cerr<<"Exception: "<<Record.Pos<<"/"<<RG.CurrentRange.first<<std::endl;
                std::__throw_invalid_argument("Using unsorted SAM file! Aborting...");
            }

            // Ignore invalid SAM records
            if(Record.ValidBit)
            {
                //Ignore empty RG 
                if(RG.size()==0);
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
                        pool.submit(IsoLasso::utils::ProcessReadGroup,std::move(RG));
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

                auto [current_start,current_end] = std::move(Record.GetRange());

                if(RG.CurrentRange.first==0)
                    RG.CurrentRange.first = current_start;

                RG.CurrentRange.second = std::max(current_end,RG.CurrentRange.second);

#ifdef DEBUG
                if(RG.validSize()<10)
                    std::cout<<"["<<RG.CurrentRange.first<<","<<RG.CurrentRange.second<<"]"<<"\n";
#endif
            }
            if(linecount%1000000==0)
            {
                const auto curtime = std::time(0);
                const auto now     = std::localtime(&curtime);
                std::cout<<"["
                         <<now->tm_year+1900<<"-"<<now->tm_mon+1<<"-"<<now->tm_mday<<" "
                         <<now->tm_hour<<":"<<now->tm_min<<":"<<now->tm_sec<<"] : "
                         <<linecount<<" records processed."<<std::endl;
            }
                
            linecount++;
        }//End of reading SamFILE

        //Last ReadGroup
        RG.RG_index = Total_RG_index;
        Valid_RG_index++;

        std::cout<<"Process last Read Group."<<std::endl;
        if(RG.ReadStart.size()>MIN_RG_SIZE)
        {
            ReadCount += RG.validSize();
#ifdef DEBUG
            std::cout<<"Processing Last ReadGroup "<<RG.RG_index<<", Range:"<<RG.ChrName<<" ["<<RG.CurrentRange.first<<","<<RG.CurrentRange.second<<"]"<<std::endl;
            std::cout<<"Number of valid reads:"<<RG.validSize()<<"\n";
#endif
            IsoLasso::utils::ProcessReadGroup(std::move(RG));
            std::cout<<"Finished last readgroup."<<std::endl;
        }
        return {Valid_RG_index,ReadCount};
    }//end of ReadSamFile
};//end of namespace IsoLaso::utils