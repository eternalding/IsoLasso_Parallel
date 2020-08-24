#ifndef ISOLASSO_FORMAT_RANGESET_HPP
#define ISOLASSO_FORMAT_RANGESET_HPP

#include <SAM_module/SAMrecord.hpp>
#include <map>

namespace IsoLasso::format
{
    class RangeSet
    {

            std::map<u_int32_t,u_int32_t>   coverage;
            bool                            Is_merged {false};

        public:

            u_int32_t
            MinDistance(const std::vector<u_int32_t>& start,const std::vector<u_int32_t>& end);

            u_int32_t
            MinDistance(const range_type& exon_range);

            u_int32_t
            MinDistance(RangeSet& SecondRange);

            inline u_int32_t
            MinDistance(const u_int32_t& pos)
            {
                std::map<u_int32_t,u_int32_t>::iterator cvg_iter {coverage.upper_bound(pos)};
                if(cvg_iter==coverage.end())
                    return abs(pos+1-(coverage.rbegin()->first));
                else if (cvg_iter->second==0)
                    return 0;
                else
                {
                    u_int32_t DistFromNext (abs(pos-cvg_iter->first));
                    if(cvg_iter!=coverage.begin())
                    {
                        cvg_iter--;
                        u_int32_t DistFromPrev (abs(pos+1-cvg_iter->first));
                        return DistFromNext>DistFromPrev?DistFromPrev:DistFromNext;
                    }
                    else   
                        return DistFromNext;
                }       
            }

            void
            Add(const std::vector<u_int32_t> start,const std::vector<u_int32_t> end);

            void
            UpdateCvg(const u_int32_t exon_start, const u_int32_t exon_end ,bool Is_merge);

            inline bool 
            IsOverlap(const u_int32_t& pos)
            {
                std::map<u_int32_t,u_int32_t>::iterator cvg_iter {coverage.upper_bound(pos)};
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

            inline std::tuple<std::vector<u_int32_t>,std::vector<u_int32_t>>
            toRange()
            {
                std::vector<u_int32_t> Start;
                std::vector<u_int32_t> End;
                
                u_int32_t prev {0};
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