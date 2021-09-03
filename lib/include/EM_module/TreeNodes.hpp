#ifndef TREENODES_HPP
#define TREENODES_HPP

#include<SAM_module/ReadGroup.hpp>

#include <iomanip>
#include <vector>
#include <queue>
#include <utility>
#include <ranges>
#include <execution>

namespace IsoLasso::Algorithm
{
    class ExonNode
    {
        private:
            std::int32_t                            Exon_index  {}  ;
            std::vector<std::pair<int32_t,int32_t>> Neighbors       ; // (Exon_index,OutDegree)
            std::int32_t                            Ancestor    {-1};
            std::int32_t                            LastVisited {-1};
        public:
            ExonNode(const std::uint32_t             Exon_index,
                     const std::vector<uint32_t>&     OutDegree,
                     const std::vector<double>&           ExpLv,
                     const std::vector<double>&          JuncLv)
            {
                auto ValidExons = std::ranges::views::iota(0,int(OutDegree.size()))
                                 |std::ranges::views::filter([&Exon_index](const uint32_t i) {return i>Exon_index;})
                                 |std::ranges::views::filter([&OutDegree](const auto i) { return OutDegree[i]>=MIN_JUNC_READ;});

                if(std::ranges::empty(ValidExons))
                    return;

                auto MaxExpLv  = std::ranges::max(ValidExons | std::ranges::views::transform([&ExpLv](auto i) -> const auto& { return ExpLv[i]; }));

                auto MaxJuncLv = std::ranges::max(ValidExons | std::ranges::views::transform([&JuncLv](auto i) -> const auto& { return JuncLv[i]; }));

                auto ValidChilds = ValidExons 
                                   |std::ranges::views::filter([&ExpLv,&JuncLv,MaxExpLv,MaxJuncLv](auto i)
                                                               {return (ExpLv[i]>=MaxExpLv*EXON_MIN_FRAC)&&(JuncLv[i]>=MaxJuncLv*JUNC_EXP_FRAC);}
                                                              )                
                                   |std::ranges::views::transform([&OutDegree](auto i){return std::pair(i,OutDegree[i]);});

                std::ranges::copy(ValidChilds,std::back_inserter(Neighbors));
                std::stable_sort(Neighbors.begin(), Neighbors.end(),
                                 [](auto i1, auto i2) {return i1.second < i2.second;});   
            }

            inline const std::vector<std::pair<int32_t,int32_t>>
            getNeighbors()
            {
                return Neighbors;
            }
            inline const std::int32_t
            getAncestor()
            {
                return Ancestor;
            }        
            inline const std::int32_t
            getLastVisited()
            {
                return LastVisited;
            }    
    };

}
#endif
