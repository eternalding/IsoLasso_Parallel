#include<SAM_module/ReadGroup.hpp>
#include<EM_module/PredExpLevel.hpp>
#include<utils/Commontype.hpp>
#include <numeric>
#include <algorithm>

namespace IsoLasso::Algorithm
{
    void
    PredExpLevel(format::ReadGroup& RG)
    {
//#ifdef DEBUG
        utils::Check_ReadGroup(RG);
//#endif

        //std::cout<<"Start Prediction for RG "<<RG.RG_index<<"-"<<RG.SubRG_index<<std::endl;
        TwoDimVec<bool>       Candidate_Isfs;
        std::vector<uint32_t> SubInsts;
        bool                  Reached_Max {true};

        //std::cout<<"Num of exons:"<<RG.ExonBoundary.size()<<std::endl;

        //Configs for prediction
        GenerateCandidateIsoform(RG,Candidate_Isfs,SubInsts);
        //std::cout<<"Stack Isoforms"<<Candidate_Isfs.size()<<std::endl;
        
        /*
        std::cout<<"SubInsts"<<std::endl;
        for(auto Insts:SubInsts)
            std::cout<<Insts<<" ";
        std::cout<<std::endl;

        std::cout<<"Candidate_Isf"<<std::endl;
        for(auto Isf:Candidate_Isfs)
        {
            for(auto idx:Isf)
                std::cout<<idx<<" ";
            std::cout<<std::endl;
        }*/

        return;
    }

    void
    GenerateCandidateIsoform(format::ReadGroup& RG,TwoDimVec<bool>& Candidate_Isfs,std::vector<uint32_t>& SubInsts)
    {
        auto NumExons {RG.ExonBoundary.size()};
        auto MaxExonExplv {0.0},MaxJuncExonExplv {0.0};
        // HyperParameters
        double exonminfrac = EXON_MIN_FRAC , 
               juncexpfrac = JUNC_EXP_FRAC ,
               intronretentionfrac = INTRON_RETENTION_FRAC;
        std::vector<bool>   IsVisited(NumExons,false),
                            IntronRetention_Left(RG.ExonBoundary.size(),false),
                            IntronRetention_Right(RG.ExonBoundary.size(),false);
        std::vector<double> Explv,
                            JuncExplv,
                            LeftExplv(RG.ExonBoundary.size(),0),
                            RightExplv(RG.ExonBoundary.size(),0);
        std::vector<uint32_t> Selected_Exon,
                              Indegree(RG.ExonBoundary.size(),0),
                              Outdegree(RG.ExonBoundary.size(),0);
        //Calculate Expression levels for each exon.
        GetExpLv(RG,Explv,JuncExplv,LeftExplv,
                 RightExplv,Indegree,Outdegree,
                 IntronRetention_Left,IntronRetention_Right,
                 MaxExonExplv,MaxJuncExonExplv);

        while(true)
        {
            bool ReachedMax {false};
            TwoDimVec<bool> CurrentPaths;

            //Choose important exons as start.
            for(auto exon_index=0;exon_index<RG.ExonBoundary.size();exon_index++)
            {
                bool Is_selected {((!(IntronRetention_Left[exon_index]||IntronRetention_Right[exon_index]) && Explv[exon_index]>=MaxExonExplv*exonminfrac )||
                                   ( (IntronRetention_Left[exon_index]||IntronRetention_Right[exon_index]) && Explv[exon_index]>=MaxExonExplv*juncexpfrac )||
                                   ( JuncExplv[exon_index] >= MaxJuncExonExplv*juncexpfrac))
                                   && RG.ExonStats[exon_index][2] <= 0.7
                                 };
                
                if(Is_selected)
                    Selected_Exon.emplace_back(exon_index); // Should start from here
            }
            std::cout<<"Selected Exons"<<std::endl;
            for(auto exon:Selected_Exon)
                std::cout<<exon<<" ";
            std::cout<<std::endl;

            while(true)
            {         
                TwoDimVec<bool>   StackIsofs;
                std::vector<bool> Path(NumExons,false);
                bool              END_FLAG {true};
                uint32_t          IsfStart,CurrentSubInst{0};    
                
                for(auto exon:Selected_Exon)
                {
                    if(!IsVisited[exon])
                    {
                        IsfStart = exon;
                        END_FLAG = false;
                        break;
                    }
                }
                //All Starts are used
                if(END_FLAG)
                    break;

                GetConnectingPaths(RG.ReadSupportMatrix,IsfStart,Path,Explv,StackIsofs,intronretentionfrac);
                
                std::cout<<"Path"<<std::endl;
                for(auto exon:Path) std::cout<<exon<<" ";
                std::cout<<std::endl;

                std::cout<<"StackIsofs"<<std::endl;
                for(auto Isf:StackIsofs) 
                {
                    for(auto exon:Isf)
                        std::cout<<exon<<" ";
                    std::cout<<std::endl;
                }
                if(Candidate_Isfs.size()+StackIsofs.size()>MAX_ISOFORM_NUM)//Exceed maximum Isoform number
                {
                    ReachedMax = true;
                    break;
                }

                auto Overlap {false};

                for(auto ExonIndex=0;ExonIndex<NumExons;ExonIndex++)//All exons
                {
                    if(Path[ExonIndex]&&IsVisited[ExonIndex])
                    {
                        Overlap = true;
                        break;
                    }
                }
                
                if(!Overlap)
                {
                    CurrentSubInst=CurrentPaths.size();
                    CurrentPaths.emplace_back(Path);
                }
                else //Overlapping
                {
                    for(auto PathIndex=0;PathIndex<CurrentPaths.size();PathIndex++)//Existing Paths
                    {
                        for(auto ExonIndex=0;ExonIndex<Path.size();ExonIndex++)//All exons
                        {
                            if(Path[ExonIndex]&&CurrentPaths[PathIndex][ExonIndex])
                            {
                                CurrentSubInst = PathIndex;
                                break;
                            }
                        }
                    }
                }

                // All StackIsofs found for current IsfStart were saved
                // to Candidate_Isfs 
                std::cout<<"StackIsofs"<<StackIsofs.size()<<std::endl;
                std::vector<uint32_t> SubInstIndex(StackIsofs.size(),CurrentSubInst);

                std::move(SubInstIndex.begin(),SubInstIndex.end(),std::back_inserter(SubInsts));
                std::move(StackIsofs.begin(),StackIsofs.end(),std::back_inserter(Candidate_Isfs));
                
                for(auto ExonIndex=0;ExonIndex<NumExons;ExonIndex++)
                    IsVisited[ExonIndex] = IsVisited[ExonIndex] || Path[ExonIndex];             
            }
            if(ReachedMax)
            {
                exonminfrac         = (exonminfrac>1)?0.99:exonminfrac*THERSHOLD_GROWTH;
                juncexpfrac         = (juncexpfrac>1)?0.99:juncexpfrac*THERSHOLD_GROWTH;
                intronretentionfrac = (intronretentionfrac>1)?0.99:intronretentionfrac*THERSHOLD_GROWTH;
                Selected_Exon.clear();
                Selected_Exon.shrink_to_fit();
                Candidate_Isfs.clear();
                Candidate_Isfs.shrink_to_fit();
                SubInsts.clear();
                SubInsts.shrink_to_fit();
                continue;
            }
            else
                break;
        }
        if(Candidate_Isfs.size()!=0)
            FilterImpossibleCandidates(Candidate_Isfs,RG);
        return;
    }

    /*
     * Get every kind of expression level
     * 
     * 
     */ 
    void 
    GetExpLv(format::ReadGroup& RG,
             std::vector<double>& Explv,
             std::vector<double>& JuncExplv,
             std::vector<double>& LeftExplv,
             std::vector<double>& RightExplv,
             std::vector<uint32_t>& Indegree,
             std::vector<uint32_t>& Outdegree,
             std::vector<bool>& IntronRetention_Left,
             std::vector<bool>& IntronRetention_Right,
             double& MaxExonExplv,
             double& MaxJuncExonExplv)
    {

        LeftExplv.assign(RG.ExonBoundary.size(),0);
        RightExplv.assign(RG.ExonBoundary.size(),0);

        for(auto row=0;row<RG.ExonBoundary.size();row++)
        {
            for(auto col=row+1;col<RG.ExonBoundary.size();col++)
            {
                if(RG.ReadSupportMatrix[row][col]>=MIN_JUNC_READ)
                {
                    Outdegree[row]++;
                    Indegree[col]++;
                }
                RightExplv[row] += RG.ReadSupportMatrix[row][col];
                LeftExplv[col]  += RG.ReadSupportMatrix[row][col];
            }
        }
        //Exon
        utils::GetExonExplv(RG,Explv);
        MaxExonExplv = *std::max_element(Explv.begin(),Explv.end());

        //Junction
        for(auto exon_index=0;exon_index<RG.ExonBoundary.size();exon_index++)
        {
            LeftExplv[exon_index]  /= (RG.ReadLen*EXON_READLEN_FRAC);
            RightExplv[exon_index] /= (RG.ReadLen*EXON_READLEN_FRAC);
            JuncExplv.emplace_back((LeftExplv[exon_index]+RightExplv[exon_index])/2);
        }
        MaxJuncExonExplv = *std::max_element(JuncExplv.begin(),JuncExplv.end());

        //IntronRetention
        for(auto exon_index=1;exon_index<RG.ExonBoundary.size()-1;exon_index++)
        {
            if(RG.ExonBoundary[exon_index].first == RG.ExonBoundary[exon_index-1].second +1)
                IntronRetention_Left.emplace_back(true);
            else
                IntronRetention_Left.emplace_back(false);
            if(RG.ExonBoundary[exon_index].second == RG.ExonBoundary[exon_index+1].first -1)
                IntronRetention_Right.emplace_back(true);
            else
                IntronRetention_Right.emplace_back(false);
        }
        return;
    }

    /*
     * Given start exon (IsfStart), get the path as an Isoform. 
     */ 

    void
    GetConnectingPaths(TwoDimVec<uint32_t>& ReadSupportMatrix,
                       uint32_t& IsfStart,
                       std::vector<bool>& Path,
                       std::vector<double>&Explv,
                       TwoDimVec<bool>& StackIsofs,
                       double& intronretentionfrac)
    {
        auto                  nExons {ReadSupportMatrix.size()};
        std::vector<uint32_t> Stack(nExons,0),Stack_Prev(nExons,0);
        std::vector<bool>     IsVisited(nExons,false);
        int32_t               Stackpt {0},NumIsofs {0};

        //Init
        Stack[0] = IsfStart;
        Stack_Prev[0] = -1;

        while(Stackpt>=0) //Keep traversing
        {
            uint32_t Current_exon {Stack[Stackpt]};
            bool     ReachedEnd   {true};
            Path[Current_exon]=true; //Visit

            if(!IsVisited[Stackpt];int32_t Prev_Stackpt=Stackpt)
            {
                IsVisited[Stackpt] = true;
                std::vector<uint32_t>& OutDegree = ReadSupportMatrix[Current_exon]; //OutDegree for current exon
                std::vector<uint32_t> ConnectedExons(OutDegree.size());
                std::iota(ConnectedExons.begin(),ConnectedExons.end(),0);
                //Sort candidate exon by OutDegree
                std::stable_sort(ConnectedExons.begin(), ConnectedExons.end(),
                                [OutDegree](auto i1, auto i2) {return OutDegree[i1] < OutDegree[i2];});

                if(OutDegree[ConnectedExons[0]]==0)//No Connected Exon
                    ConnectedExons.erase(ConnectedExons.begin(),ConnectedExons.end());
                else
                {
                    auto threshold { OutDegree[ConnectedExons[0]]<MIN_JUNC_READ
                                    ?OutDegree[ConnectedExons[0]]
                                    :MIN_JUNC_READ};
                    for(auto exon_index=0;exon_index<OutDegree.size();exon_index++)
                    {
                        if(OutDegree[ConnectedExons[exon_index]]<MIN_JUNC_READ)//OutDegree is not huge enough
                        {
                            ConnectedExons.erase(ConnectedExons.begin()+exon_index,ConnectedExons.end());
                            break;
                        }
                    }
                } 
                if(!ConnectedExons.empty()) //Candidate exons
                {
                    double MaxJuncExpLv {0.0}, MaxExpLv {0.0};
                    bool   HasCandidate {false};
                    std::vector<bool> ValidConnectedExons(ConnectedExons.size(),true);
                    for(auto index:ConnectedExons)
                    {
                        if(Explv[index]>MaxExpLv)
                            MaxExpLv = Explv[index];
                        if(ReadSupportMatrix[Current_exon][index]>MaxJuncExpLv)
                            MaxJuncExpLv = ReadSupportMatrix[Current_exon][index];
                    }
                    //Filter by Both Explv & JuncExpLv
                    for(auto index=0;index<ConnectedExons.size();index++)
                    {
                        auto candidate_exon {ConnectedExons[index]};
                        if((Explv[candidate_exon]>=MaxExpLv*EXON_MIN_FRAC) &&
                           (ReadSupportMatrix[Current_exon][candidate_exon]>=MaxJuncExpLv*JUNC_EXP_FRAC))
                            HasCandidate = true;
                        else
                            ValidConnectedExons[index] = false;
                    }
                    //If no valid candidate, use only Explv (released condition)
                    if(!HasCandidate)
                    {
                        for(auto index=0;index<ConnectedExons.size();index++)
                        {
                            auto candidate_exon {ConnectedExons[index]};
                            if(Explv[candidate_exon]>=MaxExpLv*EXON_MIN_FRAC)
                                ValidConnectedExons[index] = true;
                            else
                                ValidConnectedExons[index] = false;
                        }
                    }
                    //IntronRetention
                    auto NextExon {std::find(ConnectedExons.begin(),ConnectedExons.end(),Current_exon+1)};
                    if(NextExon!=ConnectedExons.end()) //Intron retention exist
                    {
                        if(Explv[*NextExon]<MaxExpLv*intronretentionfrac)
                            ValidConnectedExons[std::distance(ConnectedExons.begin(),NextExon)]=false;
                    }

                    //Filter exons with ValidConnectedExons = false
                    std::erase_if(ConnectedExons,[&ValidConnectedExons](auto index){return !ValidConnectedExons[index];});
                    
                    ReachedEnd = ConnectedExons.empty();

                    for(auto index:ConnectedExons)//Put all valid connected exons into stack
                    {
                        Stack.emplace_back(index);
                        Stack_Prev.emplace_back(Prev_Stackpt);
                        IsVisited.emplace_back(false);
                    }
                    Stackpt += uint32_t(ConnectedExons.size());
                }
                if(ReachedEnd)//No more candidates
                {
                    NumIsofs++;
                    int32_t  Current_Stackpt   {Stackpt}; 
                    uint32_t Current_StackExon {Current_exon};
                    std::vector<bool> Isoform(ReadSupportMatrix.size(),false);
                    while(true)
                    {
                        Isoform[Current_StackExon] = true;
                        Current_Stackpt = Stack_Prev[Current_Stackpt];
                        if(Current_Stackpt !=-1)
                            Current_StackExon = Stack[Current_Stackpt];
                        else
                            break;
                    }
                    StackIsofs.emplace_back(Isoform);
                    std::cout<<"Here"<<std::endl;
                    if(NumIsofs>=MAX_ISOFORM_NUM)
                        break;
                }
            }
            else
                Stackpt--;
        }
        return;
    }
    void
    FilterImpossibleCandidates(TwoDimVec<bool>& Candidaite_Isfs,format::ReadGroup& RG)
    {
        auto NumExons {RG.ExonBoundary.size()},NumCandidates {Candidaite_Isfs.size()};
        std::vector<bool> IsValid(NumCandidates,false);

        std::vector<std::vector<uint32_t>> SGJuncType;
        std::vector<uint32_t>              SGJuncCount;

        for(auto Typeindex=0;Typeindex<RG.SGTypes.size();Typeindex++)
        {
            if(std::accumulate(RG.SGTypes[Typeindex].begin(),RG.SGTypes[Typeindex].end(),0,[](auto i1,auto i2){return i1+(i2==true);})>1)
            {
                SGJuncType.emplace_back(RG.SGTypes[Typeindex]);
                SGJuncCount.emplace_back(RG.TypeCount[Typeindex]);
            }
        }

        for(auto Candidateindex=0;Candidateindex<NumCandidates;Candidateindex++)
        {
            if(std::accumulate(Candidaite_Isfs[Candidateindex].begin(),Candidaite_Isfs[Candidateindex].end(),0,[](auto i1,auto i2){return i1+(i2==true);})==1)
                continue;
            std::vector<bool> IsCovered(NumExons,false);
            bool InvalidFlag {false};

            for(auto Type:SGJuncType)
            {
                if(!utils::CheckCompatible(Candidaite_Isfs[Candidateindex],Type))
                    continue;
                for(auto ExonIndex=0;ExonIndex<NumExons;ExonIndex++)
                    IsCovered[ExonIndex] = IsCovered[ExonIndex] || Type[ExonIndex];
            }
            for(auto ExonIndex=0;ExonIndex<NumExons;ExonIndex++)
            {
                if(Candidaite_Isfs[Candidateindex][ExonIndex]==true && IsCovered[ExonIndex]==false)
                {
                    InvalidFlag = true;
                    break;
                }
            }
            if(InvalidFlag)
                Candidaite_Isfs[Candidateindex].assign(IsCovered.begin(),IsCovered.end());
        }


        return;
    }



}//end of IsoLasso::Algorithm

namespace IsoLasso::utils
{




}