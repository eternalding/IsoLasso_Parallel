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
#include <cmath>
#include <execution>

std::mutex IO_Mutex;
std::ofstream IsoLasso::utils::RG_STATS_FS,IsoLasso::utils::GTF_FS;
uint64_t IsoLasso::utils::TOTAL_READ_CNT;


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

        //Calculate Isoform direction
        std::vector<int64_t> IsoDir(Candidate_Isfs.size(),0);

        TwoDimVec<double> SGSupport(RG.SGTypes.size(),std::vector<double>(Candidate_Isfs.size(),0));
        GetTypeSupportAndIsfDir(RG,Candidate_Isfs,SGSupport,IsfLen,IsoDir);

        //E-M Algorithm
        EMConfig EMParameters(Candidate_Isfs.size());

#ifdef  DEBUG
        auto Start_time {std::chrono::steady_clock::now()};
#endif 
        EMAlgorithm(Candidate_Isfs,SGSupport,RG.TypeCount,EMParameters);

#ifdef  DEBUG
        IsoLasso::utils::ShowRunningTime(Start_time,"IV.","E-M Algorithm");                      
#endif

        //Write to RGFile and GTF File
        std::vector<double> ExpLv;
        //Calculate FPKM
        for(auto Isf_index=0;Isf_index<Candidate_Isfs.size();++Isf_index)
            ExpLv.push_back(EMParameters.IsoformProb[Isf_index]*double(RG.ReadCount)*10e9/(IsfLen[Isf_index]*IsoLasso::utils::TOTAL_READ_CNT));

        //Write to output
        IO_Mutex.lock();
#ifdef DEBUG
        Start_time = std::chrono::steady_clock::now();
        RG.WriteStatsToFile(IsoLasso::utils::RG_STATS_FS,Candidate_Isfs,ExpLv,IsoDir);
#endif
        RG.WritePredToGTF(IsoLasso::utils::GTF_FS,Candidate_Isfs,ExpLv,IsoDir,EMParameters.IsoformProb);

#ifdef  DEBUG
        IsoLasso::utils::ShowRunningTime(Start_time,"V.","Write to outputs");              
#endif
        IO_Mutex.unlock();
        return;
    }

    void 
    GetTypeSupportAndIsfDir(const IsoLasso::format::ReadGroup& RG,
                            const TwoDimVec<uint32_t>& Candidate_Isfs,
                            TwoDimVec<double>& SGSupport,
                            const std::vector<uint32_t>& IsfLen,
                            std::vector<int64_t>& IsfDir)
    {        
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
                    IsfDir[Isf_index] += RG.TypeDirection[Type_index];
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
        const auto nIsf = SGSupport[0].size(), nTypes = SGSupport.size();
        TwoDimVec<double> IsoEmitProb {std::move(SGSupport)}, // Initialize as SGSupport
                          Responsibility(nTypes,std::vector<double>(nIsf,0.0));

        if(EMParameters.Verbose)
        {
            std::cout<<"Init Probability:"<<std::endl;
            IsoLasso::utils::print1Dvector(EMParameters.IsoformProb);
        }

        auto LogLikelihood = 0.0, PrevLogLikelihood = LogLikelihood;

        while(Iteration<EMParameters.MaxIteration)
        {
            PrevLogLikelihood = LogLikelihood;

            //E-Step 
            LogLikelihood = EStep(IsoEmitProb,TypeCount,EMParameters,Responsibility);
            
            if(EMParameters.Verbose)
            {
                std::cout<<"IsoEmitProb"<<std::endl;
                for(int i=0;i<IsoEmitProb.size();i++)
                {
                    std::cout<<"Type "<<i<<":";
                    for(int j=0;j<IsoEmitProb[0].size();j++)
                        std::cout<<IsoEmitProb[i][j]<<" ";
                    std::cout<<std::endl;
                }
                IsoLasso::utils::print2Dvector(IsoEmitProb);
                std::cout<<"-----"<<std::endl;
                std::cout<<"Responsibility:"<<std::endl;
                for(int i=0;i<IsoEmitProb.size();i++)
                {
                    std::cout<<"Type "<<i<<":";
                    for(int j=0;j<Responsibility[0].size();j++)
                        std::cout<<Responsibility[i][j]<<" ";
                    std::cout<<std::endl;
                }
                IsoLasso::utils::print2Dvector(Responsibility);
            }
            
            MStep(IsoEmitProb,TypeCount,EMParameters,Responsibility);
            
            if(EMParameters.Verbose)
            {
                std::cout<<"["<<Iteration<<"] Loglikelihood:"<<LogLikelihood<<std::endl;
                std::cout<<"Isoform Probability:"<<std::endl;
                IsoLasso::utils::print1Dvector(EMParameters.IsoformProb);
            }

            //Stopping criteria
            if(LogLikelihood-PrevLogLikelihood<=EMParameters.MinDelta && Iteration>=EMParameters.MinIteration)
                break;
            else
            {
                IsoEmitProb = Responsibility;
                ++Iteration;
            }
        }
        return;
    }

    double
    EStep(TwoDimVec<double>& IsoEmitProb,
          const std::vector<uint32_t>& TypeCount,
          EMConfig& EMParameters,
          TwoDimVec<double>& Responsibility)
    {
        auto LogLikelihood {0.0};
        //Calculate Responsibility
        for(auto TypeIndex=0;TypeIndex<IsoEmitProb.size();++TypeIndex)
        {
            auto MarginalProb {0.0};
            for(auto IsfIndex=0;IsfIndex<IsoEmitProb[0].size();++IsfIndex)
            {
                //std::cout<<IsoEmitProb[TypeIndex][IsfIndex]<<","<<EMParameters.IsoformProb[IsfIndex]<<std::endl;
                Responsibility[TypeIndex][IsfIndex] = (IsoEmitProb[TypeIndex][IsfIndex] * EMParameters.IsoformProb[IsfIndex]);
                MarginalProb += Responsibility[TypeIndex][IsfIndex];
            }
            //If MarginalProb=0, this type cannot be generated from any Isoform!
            //Simply Set responsibility to zero
            //std::cout<<"Marginal prob:"<<MarginalProb<<std::endl;
            //std::cout<<"TypeIndex "<<TypeIndex<<std::endl;
            
            
            for_each(Responsibility[TypeIndex].begin(),
                     Responsibility[TypeIndex].end(),
                     [&MarginalProb](auto& Prob){Prob=(MarginalProb==0.0)?0.0:Prob/MarginalProb;});
            //std::cout<<"After marginalization:"<<std::endl;
            //IsoLasso::utils::print1Dvector(Responsibility[TypeIndex]);         

            LogLikelihood += (TypeCount[TypeIndex] * std::log(MarginalProb));
        }
        return LogLikelihood;
    }

    void
    MStep(TwoDimVec<double>& IsoEmitProb,
          const std::vector<uint32_t>& TypeCount,
          EMConfig& EMParameters,
          const TwoDimVec<double>& Responsibility)
    {
        //Update Isoform Probability and Isoform emission probability
        auto MarginalProb {0.0};
        for(auto IsfIndex=0;IsfIndex<IsoEmitProb[0].size();++IsfIndex)
        {
            auto IsfMarignalProb {0.0};
            for(auto TypeIndex=0;TypeIndex<IsoEmitProb.size();++TypeIndex)
            {
                IsoEmitProb[TypeIndex][IsfIndex] = TypeCount[TypeIndex]*Responsibility[TypeIndex][IsfIndex];
                IsfMarignalProb += IsoEmitProb[TypeIndex][IsfIndex];
            }
            
            for(auto TypeIndex=0;TypeIndex<IsoEmitProb.size();++TypeIndex)
                IsoEmitProb[TypeIndex][IsfIndex]= (IsfMarignalProb==0.0)?0.0:IsoEmitProb[TypeIndex][IsfIndex]/IsfMarignalProb;

            EMParameters.IsoformProb[IsfIndex] = IsfMarignalProb;
            MarginalProb += IsfMarignalProb;
        }

        for_each(EMParameters.IsoformProb.begin(),
                 EMParameters.IsoformProb.end(),
                 [&MarginalProb](auto& Prob){Prob=(MarginalProb==0.0)?0.0:Prob/MarginalProb;});

        return;        
    }

}// end of namespace IsoLasso::algorithm


