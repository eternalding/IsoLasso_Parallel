#ifndef ISOLASSO_FORMAT_READGROUP_HPP
#define ISOLASSO_FORMAT_READGROUP_HPP

#include <vector>
#include <SAM_module/SAMrecord.hpp>
#include <string>
#include <utils/Commontype.hpp>
#include <map>
#include <numeric>
#include <climits>
#include <algorithm>

namespace IsoLasso::format
{
    class ReadGroup
    {
        public:
        uint32_t                                                ReadLen             {0} ;
        std::string                                             ChrName             {""};
        TwoDimVec<uint32_t>                                     ReadStart               ;
        TwoDimVec<uint32_t>                                     ReadEnd                 ;
        std::vector<int32_t>                                    PairendTable            ;
        std::vector<bool>                                       ValidRead               ;
        std::vector<int16_t>                                    Direction               ;
        std::map<std::string,std::map<uint32_t,uint32_t>>       QNameQueryTable         ;
        std::int8_t                                             Orientation       {'.'} ;
        range_type                                              CurrentRange            ;
        std::map<uint32_t,std::vector<double>>                  CvgStats                ;
        std::vector<u_int32_t>                                  ExonCoverage            ;
        std::vector<range_type>                                 ExonBoundary            ;
        std::vector<bool>                                       ValidExons              ;
        std::vector<std::vector<uint32_t>>                      SGTypes                 ;
        std::vector<bool>                                       ValidType               ;
        std::vector<uint32_t>                                   TypeCount               ;
        std::vector<int16_t>                                    TypeDirection           ;                              
        std::vector<uint32_t>                                   Read2Type               ;

        std::uint32_t                                           RG_index            {0} ;
        std::uint32_t                                           SubRG_index         {0} ;

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
            Orientation = 0;
            ExonBoundary.clear();
            ValidExons.clear();
            CvgStats.clear();
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

        inline u_int32_t
        validSize()
        {
            return std::accumulate(ValidRead.begin(),ValidRead.end(),0,
                                    [](u_int32_t x,u_int32_t y){return x+(y==1);});
        }

        inline range_type
        getRange() const
        {
            range_type cur_range(INT_MAX,0);
            for(auto i=0;i<ReadStart.size();i++)
            {
                if(cur_range.first>ReadStart[i].front()) 
                    cur_range.first=ReadStart[i].front();
                if(cur_range.second<ReadEnd[i].back()) 
                    cur_range.second=ReadEnd[i].back();
            }
            return cur_range;
        }

        void
        SplitbyRangeSet(std::vector<ReadGroup>&,const u_int32_t& minDistance);

        void
        SplitbyDirection(std::vector<ReadGroup>&,std::vector<range_type>&);

        inline void 
        SetOrientation()
        {
            int32_t DirSum {std::accumulate(Direction.begin(),Direction.end(),0)};

            if(DirSum>0)
                Orientation = {'+'};
            else if (DirSum<0)
                Orientation = {'-'};
            return;
        }
        
        void CalculateBound(const u_int32_t MIN_JUNC_COV,const u_int32_t MIN_GAP_SPAN);
        void GetCoverage(std::map<uint32_t,int32_t>& coverage);
        void GetCvgCutPoint(const std::map<uint32_t,int32_t>&coverage,std::vector<range_type>& Cutpoint,
                            const u_int32_t threshold,const u_int32_t MIN_GAP_SPAN);
        void GetCvgStats(std::map<uint32_t,int32_t>&coverage,const std::map<uint32_t,uint32_t>& Boundary,
                        std::map<uint32_t,std::vector<double>>& CvgStat);
        void CalculateType();

        std::vector<uint32_t> 
        GetType(const std::vector<uint32_t>& SegStart,const std::vector<uint32_t>& SegEnd,
                std::vector<uint32_t>& ExonCoverage,const uint32_t& MIN_OVERLAP);
        
        void
        RemoveWeakExons(const double MIN_CVG_FRAC);
        void
        CalculateValidExons();

        inline uint32_t
        ValidRGSize()
        {
            uint32_t ValidSize {0};
            for(auto read_index=0;read_index<ReadStart.size();read_index++)
                ValidSize+=(ValidRead[read_index]==true);
            return ValidSize;
        }

        void
        WriteStatsToFile(std::ofstream& ofs);


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

    inline uint32_t
    GetOverLapping(const uint32_t& Start1,const uint32_t&End1,
                   const uint32_t& Start2,const uint32_t&End2)
    {

        uint32_t Start {std::max(Start1,Start2)}, End {std::min(End1,End2)};

        return End>Start?End-Start:0;
    }





}//end of namespace IsoLasso::utils

#endif