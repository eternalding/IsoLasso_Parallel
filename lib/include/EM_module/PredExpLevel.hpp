#ifndef PREDEXPLEVEL_HPP
#define PREDEXPLEVEL_HPP

#include<SAM_module/ReadGroup.hpp>
#include<EM_module/TreeNodes.hpp>

#include <iomanip>
#include <set>

namespace IsoLasso::Algorithm
{
    void
    PredExpLevel(format::ReadGroup& RG);

    void
    GenerateCandidateIsoform(IsoLasso::format::ReadGroup& RG,
                             TwoDimVec<uint32_t>& Candidate_Isfs,
                             std::vector<uint32_t>& SubInsts);

    void
    GetConnectingPaths(TwoDimVec<uint32_t>& ReadSupportMatrix,
                       uint32_t& IsfStart,
                       std::vector<bool>& Path,
                       std::vector<double>&Explv,
                       TwoDimVec<bool>& StackIsofs,
                       double& intronretentionfrac);

    void 
    GetExpLv(format::ReadGroup& RG,std::vector<double>& Explv,
             std::vector<double>& JuncExplv,std::vector<double>&LeftExplv,std::vector<double>&RightExplv,
             std::vector<int32_t>& Indegree,std::vector<int32_t>&Outdegree,
             std::vector<bool>& IntronRetention_Left,std::vector<bool>& IntronRetention_Right,
             double& MaxExonExplv,double& MaxJuncExonExplv);
    void
    FilterImpossibleCandidates(TwoDimVec<bool>& Candidaite_Isfs,
                               format::ReadGroup& RG);


    uint32_t
    TreeTraversal(std::vector<std::vector<uint32_t>>& TotalPaths,
                  std::vector<IsoLasso::Algorithm::ExonNode>& Tree,
                  const TwoDimVec<uint32_t>& ReadSupportMatrix,
                  const std::vector<double>& Explv,
                  const std::vector<double>& JuncExplv,
                  const uint32_t StartExon,
                  const double exon_min_frac,
                  const double junc_exp_frac ,
                  const double intron_retention_frac,
                  std::vector<bool>& ExonsCoveredByIsoform);

}//end of IsoLasso::Algorithm

namespace IsoLasso::utils
{
    inline void
    Check_ReadGroup(format::ReadGroup& RG)
    {
        std::cout<<std::setfill('=');
        std::cout<<std::right<<std::setw(40)
                 <<"Start ReadGroup "<< RG.RG_index<<" - "<<RG.SubRG_index
                 <<std::left<<std::setw(40)
                 <<" Information"<<std::endl;
        std::cout<<"ReadLength:"<<RG.ReadLen<<std::endl;
        std::cout<<"Number of exons:"<<RG.ValidExons.size()<<std::endl;
        std::cout<<"Exons:"<<std::endl;
        auto exon_index {0};
        std::cout<<"Exon"<<"\t"<<"Start"<<"\t"<<"End"<<"\t"<<"Length"<<std::endl;
        for(const auto exon:RG.ExonBoundary)
        {
            exon_index++;
            std::cout<<"E"<<exon_index<<":\t"<<exon.first<<"\t"<<exon.second<<"\t"<<exon.second-exon.first<<std::endl;
        }
        std::cout<<"SGTypes:"<<RG.SGTypes.size()<<std::endl;

        //Header for Exons
        std::cout<<"Type\t";
        for(auto idx=1;idx<=RG.ExonBoundary.size();++idx)
            std::cout<<"E"<<idx<<"\t";
        
        std::cout<<"TypeCnt"<<"\t"<<"TypeDir"<<std::endl;
        for(auto Type_idx=0;Type_idx<RG.SGTypes.size();++Type_idx)
        {
            std::cout<<Type_idx<<"\t";
            auto Exon_iter {RG.SGTypes[Type_idx].begin()};

            for(auto exon_index=0;exon_index<RG.ExonBoundary.size();++exon_index)
            {
                if(Exon_iter==RG.SGTypes[Type_idx].end())
                    std::cout<<"0\t";
                else if(*Exon_iter==exon_index)
                {
                    std::cout<<"1\t";
                    Exon_iter++;
                }
                else
                    std::cout<<"0\t";
            }            
            std::cout<<RG.TypeCount[Type_idx]<<"\t"<<RG.TypeDirection[Type_idx]<<std::endl;
        }

        std::cout<<"ReadSupportMatrix:"<<std::endl;

        for(const auto& row:RG.ReadSupportMatrix)
        {
            for(auto col:row)
                std::cout<<col<<"\t";
            std::cout<<std::endl;
        }

        std::cout<<std::right<<std::setw(40)
                 <<"End of ReadGroup " << RG.RG_index<<" - "<<RG.SubRG_index
                 <<std::left<<std::setw(40)
                 <<" Information"<<std::endl;
        return;
    }

    inline void
    GetExonWeight(IsoLasso::format::ReadGroup& RG,std::vector<double>& ExonWeight)
    {
        const double MaxExonLen {EXON_READLEN_FRAC*RG.ReadLen};
        for(const auto& Exon:RG.ExonBoundary)
        {
            auto ExonLength {Exon.second-Exon.first};
            if(ExonLength>=MaxExonLen)
                ExonWeight.push_back(1);
            else
                ExonWeight.push_back(double(ExonLength)/MaxExonLen);
        }
        return;
    }

    /*
     * Calculate the expression level of each exon 
     */
    inline void
    GetExonExplv(IsoLasso::format::ReadGroup& RG,std::vector<double>& Explv)
    {

        std::vector<double> ExonWeight;
        IsoLasso::utils::GetExonWeight(RG,ExonWeight);

        for(auto exon_index=0;exon_index<RG.ReadSupportMatrix.size();++exon_index)
        {
            auto ExonLen {RG.ExonBoundary[exon_index].second-RG.ExonBoundary[exon_index].first+1};
            auto RowSum  {std::accumulate(RG.ReadSupportMatrix[exon_index].begin(),RG.ReadSupportMatrix[exon_index].end(),0)};
            auto ColSum  {0};
            for(auto col_index=0;col_index<RG.ReadSupportMatrix[exon_index].size();col_index++)
                ColSum+=RG.ReadSupportMatrix[col_index][exon_index];

            Explv.push_back(double(RowSum+ColSum)*ExonWeight[exon_index]/ExonLen); 
        }
        return;
    }




    /*
     * Check if this Isoform can support this SGType or not
     *
     */

    inline bool
    CheckCompatible(const std::vector<bool>& Isoform,const std::vector<uint32_t>& JuncType)
    {
        const auto NumExons {Isoform.size()};
        std::vector<uint32_t> JuncIndex; 

        for(auto ExonIndex=0;ExonIndex<NumExons;++ExonIndex)//
        {
            if (Isoform[ExonIndex]==false && JuncType[ExonIndex]==true)
                return false;
            if(JuncType[ExonIndex]==true)
                JuncIndex.push_back(ExonIndex);
        }
        for(auto ExonIndex=JuncIndex.front();ExonIndex<=JuncIndex.back();ExonIndex++)
        {
            if(Isoform[ExonIndex]==true && JuncType[ExonIndex]==false)
                return false;
        }
        return true;
    }


}//end of IsoLasso::utils
#endif