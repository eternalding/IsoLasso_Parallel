#include<SAM_module/RangeSet.hpp>
#include<limits>
#include<map>
#include<cmath>
#include<algorithm>

namespace IsoLasso::format
{

    uint32_t
    RangeSet::MinDistance(const std::vector<uint32_t>& start,const std::vector<uint32_t>& end)
    {
        uint32_t current_min {std::numeric_limits<uint32_t>::max()};
        for(auto seg_index=0;seg_index<start.size();seg_index++)
        {
            uint32_t Seg_minDist {MinDistance(range_type(start[seg_index],end[seg_index]))}; //Distance for current seg to Current Rangeset
            if (Seg_minDist==0)
                return 0;
            else if (current_min> Seg_minDist)
                current_min = Seg_minDist;
        }
        return current_min;
    }

    uint32_t
    RangeSet::MinDistance(const range_type& exon_range)
    {
        auto cvg_iter {coverage.upper_bound(exon_range.first)};

        if(cvg_iter==coverage.end())//exon_range.first larger than all
            return exon_range.first+1-coverage.rbegin()->first;
        else if (cvg_iter->second==0)//Aligned on covered region
            return 0; 
        else
        {
            if(cvg_iter->first>exon_range.second)
            {
                uint32_t DistFromEnd {cvg_iter->first-exon_range.second};
                if(cvg_iter!=coverage.begin())
                {
                    cvg_iter--;
                    uint32_t DistFromStart {cvg_iter->first>exon_range.first?
                                            cvg_iter->first-exon_range.first+1:
                                            exon_range.first-cvg_iter->first+1};
                    DistFromEnd = (DistFromEnd>DistFromStart)?DistFromStart:DistFromEnd;
                }
                return DistFromEnd;
            }
            else
                return 0;
        }
    }

    uint32_t
    RangeSet::MinDistance(RangeSet& SecondRange)
    {
        uint32_t MinDist     {std::numeric_limits<uint32_t>::max()};
        uint32_t current_min {std::numeric_limits<uint32_t>::max()};

        //R1->R2
        for(auto cvg_iter=coverage.begin();cvg_iter!=coverage.end();cvg_iter++)
        {
            if(cvg_iter->second>0)//1
                current_min = SecondRange.MinDistance(cvg_iter->first);
            else
                current_min = SecondRange.MinDistance(cvg_iter->first-1);
            if (current_min==0)
                return 0;
            else
                MinDist = std::min(MinDist,current_min);
        }

        //R2->R1
        for(auto cvg_iter=SecondRange.coverage.begin();cvg_iter!=SecondRange.coverage.end();cvg_iter++)
        {
            if(cvg_iter->second>0)//1
                current_min = MinDistance(cvg_iter->first);
            else
                current_min = MinDistance(cvg_iter->first-1);
            if (current_min==0)
                return 0;
            else
                MinDist = std::min(MinDist,current_min);
        }

        return MinDist;
    }

    void
    RangeSet::Add(const std::vector<uint32_t> start,const std::vector<uint32_t> end)
    {
        for(auto i=0;i<start.size();i++)//Add each segment of the read to current RangeSet 
            UpdateCvg(start[i],end[i],false);
        merge();
        return;
    }

    void
    RangeSet::UpdateCvg(const uint32_t exon_start,const uint32_t exon_end, bool To_merge)
    {        
        bool RemoveLastElement {IsOverlap(exon_end+1)?true:false};

        coverage[exon_start]++;
        //coverage[exon_end]++;
        coverage[exon_end+1];//=0

        std::map<uint32_t,uint32_t>::iterator cvg_iter_start {coverage.find(exon_start)},
                                                cvg_iter_end {coverage.find(exon_end+1)};
        
        if(cvg_iter_start!=cvg_iter_end)
        {
            cvg_iter_start++;        
            while(cvg_iter_start!=cvg_iter_end)
                coverage.erase(cvg_iter_start++);
        }
        if(RemoveLastElement)
            coverage.erase(exon_end+1);
        return;
    }

}//end of IsoLasso::format


namespace IsoLasso::utils
{

    


}//end of IsoLasso::utils