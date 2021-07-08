#ifndef ISOLASSO_FORMAT_RANGESET_HPP
#define ISOLASSO_FORMAT_RANGESET_HPP

#include <SAM_module/SAMrecord.hpp>
#include <map>

namespace IsoLasso::format
{
    class RangeSet
    {
        public:

            std::map<uint32_t,uint32_t>   coverage;
            bool                            Is_merged {false};


            uint32_t
            MinDistance(const std::vector<uint32_t>& start,const std::vector<uint32_t>& end);

            uint32_t
            MinDistance(const range_type& exon_range);

            uint32_t
            MinDistance(RangeSet& SecondRange);

            inline uint32_t
            MinDistance(const uint32_t& pos)
            {
                std::map<uint32_t,uint32_t>::iterator cvg_iter {coverage.upper_bound(pos)};
                if(cvg_iter==coverage.end())
                    return abs(pos+1-(coverage.rbegin()->first));
                else if (cvg_iter->second==0)
                    return 0;
                else
                {
                    uint32_t DistFromNext (abs(pos-cvg_iter->first));
                    if(cvg_iter!=coverage.begin())
                    {
                        cvg_iter--;
                        uint32_t DistFromPrev (abs(pos+1-cvg_iter->first));
                        return DistFromNext>DistFromPrev?DistFromPrev:DistFromNext;
                    }
                    else   
                        return DistFromNext;
                }       
            }

            void
            Add(const std::vector<uint32_t> start,const std::vector<uint32_t> end);

            void
            UpdateCvg(const uint32_t exon_start, const uint32_t exon_end ,bool Is_merge);

            inline bool 
            IsOverlap(const uint32_t& pos)
            {
                std::map<uint32_t,uint32_t>::iterator cvg_iter {coverage.upper_bound(pos)};
                if (cvg_iter==coverage.end())
                    return false;
                else if(cvg_iter->second==0)
                    return true;
                else
                    return false;
            }

            inline void
            merge()
            {
                bool Is_covered {false};
                for(auto cvg_iter=coverage.begin();cvg_iter!=coverage.end();)
                {
                    if(Is_covered!=(cvg_iter->second>0))
                    {
                        Is_covered = (cvg_iter->second>0);
                        cvg_iter->second = (cvg_iter->second>0?1:0);
                        cvg_iter++;
                    }
                    else
                        coverage.erase(cvg_iter++);
                }
                Is_merged = true;
                return;
            }

            inline std::tuple<std::vector<uint32_t>,std::vector<uint32_t>>
            toRange()
            {
                std::vector<uint32_t> Start;
                std::vector<uint32_t> End;
                
                uint32_t prev {0};
                for(auto cvg_iter=coverage.begin();cvg_iter!=coverage.end();cvg_iter++)
                {
                    if(cvg_iter->second>0)
                        prev = cvg_iter->first;
                    else
                    {
                        Start.emplace_back(prev);
                        End.emplace_back(cvg_iter->first-1);
                    }
                }
                return std::make_tuple(Start,End);
            }

    };

}//end of IsoLasso::format

namespace IsoLasso::utils
{





}//end of IsoLasso::utils


#endif