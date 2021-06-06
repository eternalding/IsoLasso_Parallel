#include <SAM_module/SAMrecord.hpp>
#include <SAM_module/SAMprocess.hpp>
#include <SAM_module/ReadGroup.hpp>
#include <utils/Commontype.hpp>
#include <utils/ArgParser.hpp>
#include <utils/Auxiliario.hpp>
#include <iostream>
#include <fstream>
#include <numeric>


namespace IsoLasso::utils
{
    std::tuple<uint32_t,uint32_t>
    ReadSamFile(const IsoLasso::format::Args& arguments)
    {
        std::ifstream fin(arguments.SAMFILE,std::ios::in);
        if (!fin.is_open())
            throw std::invalid_argument("Failed to open target SAM file!");

        IsoLasso::format::Header_record HRecord;
        IsoLasso::format::Sam_record Record;
        IsoLasso::format::ReadGroup RG;
        uint32_t Total_RG_index  {0};
        uint32_t Valid_RG_index  {0};
        uint32_t ReadCount       {0};
        std::vector<uint32_t> rg_size;

        //Header is now unavailable
        while(fin>>HRecord);
        
        //Main workflow
        while(fin>>Record)
        {
            //SAM file should be sorted
            if(Record.Pos<RG.CurrentRange.first)
                std::__throw_invalid_argument("Using unsorted SAM file! Aborting...");

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
                        ProcessReadGroup(RG);
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
        }

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
        return {Valid_RG_index,ReadCount};
    }//end of ReadSamFile
};//end of namespace IsoLaso::utils