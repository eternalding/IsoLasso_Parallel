#include<SAM_module/ReadGroup.hpp>
#include<EM_module/PredExpLevel.hpp>
#include<EM_module/TreeNodes.hpp>
#include<EM_module/EMalgorithm.hpp>
#include<utils/Commontype.hpp>
#include<numeric>
#include<algorithm>
#include<stack>
#include<set>
#include<utils/Auxiliario.hpp>
#include<ranges>
#include <mutex>

std::mutex IO_Mutex;
std::ofstream IsoLasso::utils::RG_STATS_FS,IsoLasso::utils::GTF_FS;

namespace IsoLasso::Algorithm
{

    void
    EM_Process(IsoLasso::format::ReadGroup& RG,
               TwoDimVec<uint32_t>& Candidate_Isfs,
               std::vector<uint32_t>& SubInsts)
    {
        //Calculate isoform length
        std::vector<uint32_t> IsfLen;
        IsoLasso::Algorithm::GetIsoformLen(RG,Candidate_Isfs,IsfLen);

        //Calculate Isoform direction(TODO)
        std::vector<int16_t> IsoDir;
        std::vector<double> Priors;

        TwoDimVec<double> SGSupport(RG.SGTypes.size(),std::vector<double>(Candidate_Isfs.size(),0));
        GetTypeSupport(RG,Candidate_Isfs,SGSupport,IsfLen);

        /*
        std::cout<<"Assembled Isoforms:"<<std::endl;
        IsoLasso::utils::print2Dvector(Candidate_Isfs);

        std::cout<<"TypeSupport:"<<std::endl;
        IsoLasso::utils::print2Dvector(SGSupport);
        */
        
        //E-M Algorithm
        EMConfig EMParameters(Candidate_Isfs.size());
        EMAlgorithm(Candidate_Isfs,SGSupport,RG.TypeCount,EMParameters);

        //Write to RGFile and GTF File
        std::vector<double> ExpLv;

        //Calculate FPKM
        for(auto Isf_index=0;Isf_index<Candidate_Isfs.size();++Isf_index)
            ExpLv.emplace_back(EMParameters.IsoformProb[Isf_index]*double(RG.ReadCount)/IsfLen[Isf_index]);

        IO_Mutex.lock();
        RG.WriteStatsToFile(IsoLasso::utils::RG_STATS_FS,Candidate_Isfs,ExpLv);
        RG.WritePredToGTF(IsoLasso::utils::GTF_FS,Candidate_Isfs,ExpLv);
        IO_Mutex.unlock();
        return;
    }
    void 
    GetTypeSupport(const IsoLasso::format::ReadGroup& RG,
                   const TwoDimVec<uint32_t>& Candidate_Isfs,
                   TwoDimVec<double>& SGSupport,
                   const std::vector<uint32_t>& IsfLen)
    {
        const auto& ExonBoundary = RG.ExonBoundary;
        for(auto Isf_index=0;Isf_index<Candidate_Isfs.size();Isf_index++)
        {
            for(auto Type_index=0;Type_index<RG.SGTypes.size();Type_index++)
            {
                if(CheckCompatible(Candidate_Isfs[Isf_index],RG.SGTypes[Type_index]))
                {  
                    auto TypeLen {0};
                    for(const auto ExonIndex:RG.SGTypes[Type_index])
                        TypeLen += (RG.ExonBoundary[ExonIndex].second-RG.ExonBoundary[ExonIndex].first+1);

                    SGSupport[Type_index][Isf_index] = TypeLen / double(IsfLen[Isf_index]);
                }
            }
        }
        return;
    }

    void
    EMAlgorithm(const TwoDimVec<uint32_t>& Candidate_Isfs,
                const TwoDimVec<double>& SGSupport,
                const std::vector<uint32_t>& TypeCount,
                EMConfig& EMParameters)
    {
        auto Iteration {0};
        auto delta   {0.0};
        TwoDimVec<double> IsoEmitProb = SGSupport,Responsibility(SGSupport.size(),std::vector<double>(SGSupport[0].size(),0.0)); // Initialize as SGSupport

        const auto TotalReadCnt {std::accumulate(TypeCount.begin(),TypeCount.end(),0)};

        std::vector<double> Prev_IsoformProb;

        while(Iteration<EMParameters.MaxIteration)
        {
            Prev_IsoformProb = EMParameters.IsoformProb;

            //E-Step 
            double JointProbability = EStep(IsoEmitProb,TypeCount,EMParameters,Responsibility);
            /*if(EMParameters.Verbose)
            {
                std::cout<<"IsoEmitProb"<<std::endl;
                IsoLasso::utils::print2Dvector(IsoEmitProb);

                std::cout<<"Responsibility:"<<std::endl;
                IsoLasso::utils::print2Dvector(Responsibility);
            }*/
            
            MStep(IsoEmitProb,TypeCount,EMParameters,Responsibility,TotalReadCnt);
            
            if(EMParameters.Verbose)
            {
                std::cout<<"["<<Iteration<<"]Isoform Probability:"<<std::endl;
                for(auto prob:EMParameters.IsoformProb)
                    std::cout<<prob<<" ";
                std::cout<<std::endl;
            }

            //Stopping criteria
            for(auto Isf_idx=0;Isf_idx<EMParameters.IsoformProb.size();++Isf_idx)
                delta+=(std::abs(EMParameters.IsoformProb[Isf_idx]-Prev_IsoformProb[Isf_idx]));
            if(delta<=EMParameters.MinDelta && Iteration>=EMParameters.MinIteration)
                break;
            else
                delta = 0.0;
            ++Iteration;
        }
    }

    double
    EStep(const TwoDimVec<double>& IsoEmitProb,
          const std::vector<uint32_t>& TypeCount,
          EMConfig& EMParameters,
          TwoDimVec<double>& Responsibility)
    {
        //Calculate Responsibility
        for(auto TypeIndex=0;TypeIndex<IsoEmitProb.size();++TypeIndex)
        {
            auto MarginalProb {0.0};
            for(auto IsfIndex=0;IsfIndex<IsoEmitProb[0].size();++IsfIndex)
            {
                Responsibility[TypeIndex][IsfIndex] = (IsoEmitProb[TypeIndex][IsfIndex]+TypeCount[TypeIndex]* EMParameters.bias)*EMParameters.IsoformProb[IsfIndex];
                MarginalProb += Responsibility[TypeIndex][IsfIndex];
            }
            //If MarginalProb=0, this type cannot be generated from any Isoform!
            //Simply Set responsibility to zero
            for_each(Responsibility[TypeIndex].begin(),Responsibility[TypeIndex].end(),[&MarginalProb](auto& Prob){Prob=(MarginalProb==0.0)?0.0:Prob/MarginalProb;});
        }
        return 0.0;
    }

    void
    MStep(TwoDimVec<double>& IsoEmitProb,
          const std::vector<uint32_t>& TypeCount,
          EMConfig& EMParameters,
          const TwoDimVec<double>& Responsibility,
          const uint32_t TotalReadCnt)
    {
        //Update Isoform Probability and Isoform emission probability
        for(auto IsfIndex=0;IsfIndex<IsoEmitProb[0].size();++IsfIndex)
        {
            auto MarginalProb {0};
            for(auto TypeIndex=0;TypeIndex<IsoEmitProb.size();++TypeIndex)
            {
                IsoEmitProb[TypeIndex][IsfIndex] = TypeCount[TypeIndex]*Responsibility[TypeIndex][IsfIndex];
                MarginalProb += IsoEmitProb[TypeIndex][IsfIndex];
            }
            //std::cout<<"MarginalProb:"<<MarginalProb<<std::endl;

            //If MarginalProb=0, no type can be generated from curring Isoform
            //Simply Set responsibility to zero            
            EMParameters.IsoformProb[IsfIndex] = MarginalProb/double(TotalReadCnt);

            for(auto TypeIndex=0;TypeIndex<IsoEmitProb.size();++TypeIndex)
                IsoEmitProb[TypeIndex][IsfIndex]= (MarginalProb==0.0)?0.0:IsoEmitProb[TypeIndex][IsfIndex]/MarginalProb;
        }
        return;        
    }

}// end of namespace IsoLasso::algorithm


