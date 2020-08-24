#include<SAM_module/RangeSet.hpp>
#include<limits>
#include<map>
#include<cmath>


namespace IsoLasso::format
{

    u_int32_t
    RangeSet::MinDistance(const std::vector<u_int32_t>& start,const std::vector<u_int32_t>& end)
    {
        u_int32_t current_min {std::numeric_limits<u_int32_t>::max()};
        for(auto exon_index=0;exon_index<start.size();exon_index++)
        {
            u_long Exon_minDist {MinDistance(range_type(start[exon_index],end[exon_index]))};
            if (Exon_minDist==0)
                return 0;
            else if (current_min> Exon_minDist)
                current_min = Exon_minDist;
        }
        return current_min;
    }

    u_int32_t
    RangeSet::MinDistance(const range_type& exon_range)
    {
        std::map<u_int32_t,u_int32_t>::iterator cvg_iter {coverage.upper_bound(exon_range.first)};

        if(cvg_iter==coverage.end())
            return exon_range.first+1-coverage.rbegin()->first;
        else if (cvg_iter->second==0)
            return 0; 
        else
        {
            if(cvg_iter->first>exon_range.second)
            {
                u_int32_t DistFromEnd {int32_t(cvg_iter->first)-int32_t(exon_range.second)};
                if(cvg_iter!=coverage.begin())
                {
                    cvg_iter--;
                    u_int32_t DistFromStart {std::abs(int32_t(exon_range.first)+1-int32_t(cvg_iter->first))};
                    DistFromEnd = (DistFromEnd>DistFromStart)?DistFromStart:DistFromEnd;
                }
                return DistFromEnd;
            }
            else
                return 0;
        }
    }

    u_int32_t
    RangeSet::MinDistance(RangeSet& SecondRange)
    {
        u_int32_t MinDist {std::numeric_limits<u_int32_t>::max()};
        u_int32_t current_min {std::numeric_limits<u_int32_t>::max()};

        //R1->R2
        for(auto cvg_iter=coverage.begin();cvg_iter!=coverage.end();cvg_iter++)
        {
            if(cvg_iter->second>0)
                current_min = SecondRange.MinDistance(cvg_iter->first);
            else
                current_min = SecondRange.MinDistance(cvg_iter->first-1);
            if (current_min==0)
                return 0;
            else
                MinDist = MinDist>current_min?current_min:MinDist;
        }

        //R2->R1
        for(auto cvg_iter=SecondRange.coverage.begin();cvg_iter!=SecondRange.coverage.end();cvg_iter++)
        {
            if(cvg_iter->second>0)
                current_min = MinDistance(cvg_iter->first);
            else
                current_min = MinDistance(cvg_iter->first-1);
            if (current_min==0)
                return 0;
            else
                MinDist = MinDist>current_min?current_min:MinDist;
        }

        return MinDist;
    }

    void
    RangeSet::Add(const std::vector<u_int32_t> start,const std::vector<u_int32_t> end)
    {
        for(auto i=0;i<start.size();i++)
            UpdateCvg(start[i],end[i],false);
        merge();
        return;
    }

    void
    RangeSet::UpdateCvg(const u_int32_t exon_start,const u_int32_t exon_end, bool To_merge)
    {        
        bool RemoveLastElement {IsOverlap(exon_end+1)?true:false};

        coverage[exon_start]++;
        coverage[exon_end]++;
        coverage[exon_end+1];

        std::map<u_int32_t,u_int32_t>::iterator cvg_iter_start {coverage.find(exon_start)},\
                                                cvg_iter_end   {coverage.find(exon_end)};
        
        if(cvg_iter_start!=cvg_iter_end)
        {
            cvg_iter_start++;        
            while(cvg_iter_start!=cvg_iter_end)
                coverage.erase(cvg_iter_start++);
        }
        if(RemoveLastElement)
            coverage.erase(exon_end+1);
        Is_merged=false;
        return;
    }

}//end of IsoLasso::format


namespace IsoLasso::utils
{

    


}//end of IsoLasso::utils