#ifndef ISOLASSO_FORMAT_READGROUP_HPP
#define ISOLASSO_FORMAT_READGROUP_HPP

#include <vector>
#include <SAM_module/SAMrecord.hpp>
#include <string>
#include <utils/Commontype.hpp>
#include <map>
#include <numeric>

namespace IsoLasso::format
{
    class ReadGroup
    {
        public:
        uint32_t                                                ReadLen             {0} ;
        std::string                                             ChrName             {""};
        TwoDimVec<uint32_t>                                     ReadStart               ;
        TwoDimVec<uint32_t>                                     ReadEnd                 ;
        std::vector<int32_t>                                   PairendTable            ;
        std::vector<bool>                                       ValidRead               ;
        std::vector<uint32_t>                                   Direction               ;
        std::map<std::string,std::map<uint32_t,uint32_t>>       QNameQueryTable         ;

        inline void 
        reset()
        {
            ReadStart.resize(0);
            ReadEnd.resize(0);
            PairendTable.resize(0);
            ValidRead.resize(0);
            Direction.resize(0);
            ReadStart.shrink_to_fit();
            ReadEnd.shrink_to_fit();
            PairendTable.shrink_to_fit();
            ValidRead.shrink_to_fit();
            Direction.shrink_to_fit();
            ReadLen = 0;
            ChrName = "";
            QNameQueryTable.clear();
            return;
        }

        void 
        AddRecord(const IsoLasso::format::Sam_record& record,range_type& current_range);

        inline void
        AddPair(const std::vector<u_int32_t>& Start, const std::vector<u_int32_t>& End, const u_int16_t direction,\
                const std::vector<u_int32_t>& PE_Start, const std::vector<u_int32_t>& PE_End, const u_int16_t PE_direction)
        {
            AddWithoutPair(Start,End,direction);
            AddWithoutPair(PE_Start,PE_End,PE_direction);
            PairendTable[PairendTable.size()-1] = PairendTable.size()-2;
            PairendTable[PairendTable.size()-2] = PairendTable.size()-1;
            return;
        }

        void
        RemoveLongSpanReads(const uint32_t& MAX_EXON_SPAN,const uint32_t& MAX_JUNCTION_COVERAGE);        

        inline void
        AddWithoutPair(const IsoLasso::format::Sam_record& Record)
        {
            ReadStart.emplace_back(Record.SegmentStart);
            ReadEnd.emplace_back(Record.SegmentEnd);
            Direction.emplace_back(Record.SpliceDir);
            ValidRead.emplace_back(true);
            PairendTable.emplace_back(-1);
            return;
        }

        inline void
        AddWithoutPair(const std::vector<u_int32_t>& Start, const std::vector<u_int32_t>& End, const u_int16_t direction)
        {
            ReadStart.emplace_back(Start);
            ReadEnd.emplace_back(End);
            Direction.emplace_back(direction);
            ValidRead.emplace_back(true);
            PairendTable.emplace_back(-1);
            return;
        }



        void
        SplitbyRangeSet(std::vector<IsoLasso::format::ReadGroup>& SubRGs,const uint32_t MIN_GAP_SPAN);

        inline u_int32_t
        validSize()
        {
            return std::accumulate(ValidRead.begin(),ValidRead.end(),0,
                                    [](u_int32_t x,u_int32_t y){return x+(y==1);});
        }



    };


}//end of namespace IsoLasso::format

namespace IsoLasso::utils
{
    void 
    ProcessReadGroup(IsoLasso::format::ReadGroup& RG,const range_type& current_range);

    void
    GenerateInstance(IsoLasso::format::ReadGroup& RG,const range_type& current_range);

    inline bool
    ExceedThreshold(std::map<u_int32_t,u_int32_t>& ReadCov,const u_int32_t prev_end,const u_int32_t curr_start,uint32_t threshold)
    {   
        std::map<u_int32_t,u_int32_t>::const_iterator iter = ReadCov.lower_bound(prev_end);

        uint32_t read_count {0};

        while(iter!= ReadCov.end() && iter->first<curr_start)
        {
            read_count+=iter->second;
            if(read_count>threshold)
                return true;
            else
                iter++;
        }

        return false;
    }




}//end of namespace IsoLasso::utils

#endif