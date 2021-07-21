#ifndef EMALGORITHM_HPP
#define EMALGORITHM_HPP

#include<SAM_module/ReadGroup.hpp>
#include<utils/Commontype.hpp>
#include<algorithm>
#include<ranges>
#include<vector>
#include<string>

namespace IsoLasso::Algorithm
{

    struct EMConfig
    {
        //HyperParameters 
        static constexpr   double           CutRate      {1e-6};
        static constexpr   double           MinPICut     {1e-2};
        static constexpr   double           AlphaLow     {0.0} ;
        static constexpr   double           AlphaHigh    {5.0} ;
        static constexpr   std::string_view InitMethod   {"UNIFORM"}; //RANDOM,UNIFORM
        static constexpr   uint32_t         MaxIteration {500};
        static constexpr   uint32_t         MinIteration {10};
        static constexpr   bool             Verbose      {false};
        static constexpr   double           MinDelta     {1e-4};
        static constexpr   double           bias         {0.0};
        
        //Learnable parameters
                    double  Tao        {0.0} ;
        std::vector<double> IsoformProb      ;
        

        EMConfig(const auto NumofIsoform)
        {
            if(InitMethod=="RANDOM")
            {
                IsoformProb = std::vector<double>(NumofIsoform);
                std::generate_n(IsoformProb.begin(),NumofIsoform,[](){return rand()/(double)RAND_MAX;});
                auto MarginalProb {std::accumulate(IsoformProb.begin(),IsoformProb.end(),0.0)};
                for_each(IsoformProb.begin(),IsoformProb.end(),
                        [&MarginalProb](auto& Prob){Prob/=MarginalProb;});
            }
            else if(InitMethod=="UNIFORM")     
                IsoformProb = std::vector<double>(NumofIsoform,1/double(NumofIsoform));
        }

    };

    void
    EM_Process(IsoLasso::format::ReadGroup& RG,
               TwoDimVec<uint32_t>& Candidate_Isfs,
               std::vector<uint32_t>& SubInsts);

    inline void
    GetIsoformLen(IsoLasso::format::ReadGroup& RG,
                  TwoDimVec<uint32_t>& Candidate_Isfs,
                  std::vector<uint32_t>& IsoformLen)
    {
        for(const auto& Isf:Candidate_Isfs)
        {
            auto IsfLen {0};
            for(const auto exon_index:Isf)
                IsfLen += (RG.ExonBoundary[exon_index].second-RG.ExonBoundary[exon_index].first+1);
            IsfLen = IsfLen - RG.ReadLen + 1;
            IsoformLen.push_back(IsfLen);
        }
        return;
    }

    inline bool
    CheckCompatible(const std::vector<uint32_t>& candIsf,const std::vector<uint32_t>&SGType)
    {
        if(SGType.empty())
            return false;
        else if(!std::ranges::includes(candIsf,SGType))
            return false;
        return true;
    }

    inline void
    GetTypeSupportAndIsfDir(const IsoLasso::format::ReadGroup& RG,
                            const TwoDimVec<uint32_t>& Candidate_Isfs,
                            TwoDimVec<double>& SGSupport,
                            const std::vector<uint32_t>& IsfLen,
                            std::vector<int64_t>& IsfDir);

    void
    EMAlgorithm(const TwoDimVec<uint32_t>& Candidate_Isfs,
                const TwoDimVec<double>& SGSupport,
                const std::vector<uint32_t>& TypeCount,
                EMConfig& EMParameters);
    double
    EStep(const TwoDimVec<double>& IsoEmitProb,
          const std::vector<uint32_t>& TypeCount,
          EMConfig& EMParameters,
          TwoDimVec<double>& Responsibility);

    void
    MStep(TwoDimVec<double>& IsoEmitProb,
          const std::vector<uint32_t>& TypeCount,
          EMConfig& EMParameters,
          const TwoDimVec<double>& Responsibility,
          const uint32_t TotalReadCnt);

}



#endif
