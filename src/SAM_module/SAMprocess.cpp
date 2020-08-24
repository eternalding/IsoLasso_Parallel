#include <SAM_module/SAMrecord.hpp>
#include <SAM_module/SAMprocess.hpp>
#include <SAM_module/ReadGroup.hpp>
#include <utils/Commontype.hpp>
#include <utils/ArgParser.hpp>
#include <iostream>
#include <fstream>
#include <numeric>


namespace IsoLasso::utils
{
    std::tuple<u_int32_t,u_int32_t>
    ReadSamFile(const IsoLasso::format::Args& arguments)
    {
        std::ifstream fin(arguments.SAMFILE,std::ios::in);
        if (!fin.is_open())
            throw std::invalid_argument("Failed to open target SAM file!");

        IsoLasso::format::Header_record HRecord;
        IsoLasso::format::Sam_record Record;
        IsoLasso::format::ReadGroup RG;
        range_type CurrentRange(0,0);


        u_int32_t RG_index {0};
        u_int32_t ReadCount {0};

        //Header is now unavailable
        while(fin>>HRecord);
        
        std::vector<u_int32_t> rg_size;

        //Main workflow
        while(fin>>Record)
        {
            if(Record.Pos<CurrentRange.first)
                std::__throw_invalid_argument("Using unsorted SAM file! Aborting...");

            if(Record.ValidBit)
            {
                //New ReadGroup
                if(Record.RName!=RG.ChrName || Record.Pos > CurrentRange.second+MIN_GAP_SPAN)
                {

                    RG_index ++;
                    rg_size.emplace_back(RG.ReadStart.size());

                    ReadCount+=RG.validSize();                     
                    if(RG.ReadStart.size()>MIN_RG_SIZE)
                    {
                        ProcessReadGroup(RG,CurrentRange);
                        RG.reset();
                    }
                    else if(RG.ReadStart.size()>=0) // RG is not large enough
                        RG.reset();                    
                }

                //Initialize ReadGroup
                if(RG.ChrName=="")
                {
                    RG.ReadLen=GetEfficientLen(Record);
                    RG.ChrName=Record.RName;
                }

                //Add current record to current ReadGroup
                RG.AddRecord(Record,CurrentRange);
                u_int32_t current_end {Record.GetRange().second};

                CurrentRange.second = current_end>CurrentRange.second?current_end:CurrentRange.second;
            }
        }

        RG_index ++;
        ReadCount+=RG.validSize(); 

        return {RG_index,ReadCount};
    }//end of ReadSamFile


    

};//end of namespace IsoLaso::utils