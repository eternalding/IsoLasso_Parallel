#ifndef ISOLASSO_FORMAT_READGROUP_HPP
#define ISOLASSO_FORMAT_READGROUP_HPP

#include <vector>
#include <SAM_module/SAMrecord.hpp>
#include <string>
#include <utils/Auxiliario.hpp>
#include <utils/Commontype.hpp>
#include <map>
#include <numeric>
#include <climits>
#include <algorithm>
#include <limits>

namespace IsoLasso::format
{
    class ReadGroup
    {
        public:
        uint32_t                                                ReadLen           {0}                   ;
        std::string                                             ChrName           {""}                  ;
        TwoDimVec<uint32_t>                                     ReadStart                               ;
        TwoDimVec<uint32_t>                                     ReadEnd                                 ;
        std::vector<int32_t>                                    PairendTable                            ;
        std::vector<bool>                                       ValidRead                               ;
        std::vector<int64_t>                                    Direction                               ;
        std::map<std::string,std::map<uint32_t,uint32_t>>       QNameQueryTable                         ;
        std::int8_t                                             Orientation       {'.'}                 ;
        range_type                                              CurrentRange      {range_type(0,0)}     ;
        std::map<uint32_t,std::vector<double>>                  CvgStats                                ;
        std::vector<uint32_t>                                   ExonCoverage                            ;
        std::vector<range_type>                                 ExonBoundary                            ;
        std::vector<bool>                                       ValidExons                              ;
        std::vector<std::vector<uint32_t>>                      SGTypes                                 ;
        std::vector<bool>                                       ValidType                               ;
        std::vector<uint32_t>                                   TypeCount                               ;
        std::vector<int64_t>                                    TypeDirection                           ;                              
        std::vector<uint32_t>                                   Read2Type                               ;
        std::map<uint32_t,uint32_t>                             ReadLen_Count                           ;
        std::vector<std::vector<double>>                        ExonStats                               ;
        uint32_t                                                ReadCount         {0}                   ;
        uint32_t                                                PEReadCount       {0}                   ;

        std::vector<std::vector<uint32_t>>                      PETypes                                 ;

        std::uint32_t                                           RG_index          {0}                   ;
        std::uint32_t                                           SubRG_index       {0}                   ;
        std::vector<std::vector<uint32_t>>                      ReadSupportMatrix                       ;
        std::uint32_t                                           PE_READ_DISTANCE  {PE_READ_DISTANCE}    ;
        std::uint32_t                                           PE_READ_STD       {PE_READ_STD}         ;

        inline void 
        reset()
        {
            ReadStart.resize(0);
            ReadEnd.resize(0);
            PairendTable.resize(0);
            ValidRead.resize(0);
            Direction.resize(0);
            TypeDirection.resize(0);
            TypeDirection.shrink_to_fit();
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
            CurrentRange=range_type(0,0);
            ReadLen_Count.clear();
            return;
        }

        void 
        AddRecord(const IsoLasso::format::Sam_record& record);

        inline void
        AddPair(const std::vector<uint32_t>& Start, const std::vector<uint32_t>& End, const int64_t direction,\
                const std::vector<uint32_t>& PE_Start, const std::vector<uint32_t>& PE_End, const int64_t PE_direction)
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
        AddWithoutPair(const std::vector<uint32_t>& Start, 
                       const std::vector<uint32_t>& End, 
                       const int64_t direction)
        {
            ReadStart.emplace_back(Start);
            ReadEnd.emplace_back(End);
            Direction.emplace_back(direction);
            ValidRead.emplace_back(true);
            PairendTable.emplace_back(-1);
            return;
        }

        inline uint32_t
        validSize()
        {
            return std::accumulate(ValidRead.begin(),
                                   ValidRead.end(),
                                   0,
                                   [](uint32_t x,uint32_t y){return x+(y==1);});
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
        SplitbyRangeSet(std::vector<ReadGroup>&,const uint32_t& minDistance);

        void
        SplitbyDirection(std::vector<ReadGroup>&,std::vector<range_type>&);

        inline void 
        SetOrientation()
        {
            int64_t DirSum {std::accumulate(Direction.begin(),Direction.end(),0)};

            if(DirSum>0)
                Orientation = {'+'};
            else if (DirSum<0)
                Orientation = {'-'};
            return;
        }
        
        void 
        CalculateBound(const uint32_t MIN_JUNC_COV,const uint32_t MIN_GAP_SPAN);
        void 
        GetCoverage(std::map<uint32_t,int32_t>& coverage);
        void 
        GetCvgCutPoint(const std::map<uint32_t,int32_t>&coverage,std::vector<range_type>& Cutpoint,
                       const uint32_t threshold,const uint32_t MIN_GAP_SPAN);
        void 
        GetCvgStats(std::map<uint32_t,int32_t>&coverage,const std::map<uint32_t,uint32_t>& Boundary,
                    std::map<uint32_t,std::vector<double>>& CvgStat);
        void 
        CalculateType();

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
        WriteStatsToFile(std::ofstream& ofs,
                         const TwoDimVec<uint32_t>& CandidateIsf,
                         const std::vector<double>& ExpLv,
                         const std::vector<int64_t>& IsfDir);

        void
        WritePredToGTF(std::ofstream& ofs,
                       const TwoDimVec<uint32_t>& CandidateIsf,
                       const std::vector<double>& ExpLv,
                       const std::vector<int64_t>& IsfDir);

        void
        PostProcess();

        inline void
        fillRSMatrix()
        {
            ReadSupportMatrix.assign(ExonBoundary.size(),std::vector<uint32_t>(ExonBoundary.size(),0));
            for(auto Type_index=0;Type_index<SGTypes.size();Type_index++)
            {
                if(SGTypes[Type_index].size()==1) //Single exon type
                    ReadSupportMatrix[SGTypes[Type_index][0]][SGTypes[Type_index][0]] = TypeCount[Type_index];
                else if (SGTypes[Type_index].size()>1) //Multiple exons type
                {
                    for(auto exon_index=0;exon_index<SGTypes[Type_index].size()-1;exon_index++)
                        ReadSupportMatrix[SGTypes[Type_index][exon_index]][SGTypes[Type_index][exon_index+1]] += TypeCount[Type_index];
                }
            }
            return;
        }

        inline auto
        size()
        {
            return ValidRead.size();
        }





    };


}//end of namespace IsoLasso::format

namespace IsoLasso::utils
{
    void
    ProcessReadGroup(IsoLasso::format::ReadGroup RG);

    void
    GenerateInstance(IsoLasso::format::ReadGroup& RG,const range_type& current_range);

    //Check how many reads were aligned to the junction between prev_end & curr_start
    inline bool
    ExceedThreshold(std::map<uint32_t,uint32_t>& ReadCov,const uint32_t prev_end,
                    const uint32_t curr_start,uint32_t threshold)
    {   
        std::map<uint32_t,uint32_t>::const_iterator iter = ReadCov.lower_bound(prev_end);

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

        return End>(Start-1)?End-Start+1:0;
    }

    //Get shortest distance to Target in boundary 
    inline auto
    ShortestDist(const auto& Boundary,const auto Target)
    {
        auto lower_iter = Boundary.lower_bound(Target); 
        if( lower_iter==Boundary.end())
            return std::numeric_limits<uint32_t>::max();
        if (lower_iter == Boundary.begin())
            return lower_iter->first - Target;
        const auto Upper_iter = std::prev(lower_iter); 
        if ((Target-Upper_iter->first) >= (lower_iter->first-Target))
            return lower_iter->first-Target;   
        return Target-Upper_iter->first;
    }
}//end of namespace IsoLasso::utils

#endif