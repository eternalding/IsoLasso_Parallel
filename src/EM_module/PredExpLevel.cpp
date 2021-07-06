#include<SAM_module/ReadGroup.hpp>
#include<EM_module/PredExpLevel.hpp>
#include<EM_module/TreeNodes.hpp>
#include<utils/Commontype.hpp>
#include<numeric>
#include<algorithm>
#include<stack>
#include<set>
#include<utils/Auxiliario.hpp>
#include<ranges>


namespace IsoLasso::Algorithm
{
    void
    PredExpLevel(format::ReadGroup& RG)
    {
#ifdef DEBUG
        utils::Check_ReadGroup(RG);
#endif

        TwoDimVec<bool>       Candidate_Isfs;
        std::vector<uint32_t> SubInsts;
        bool                  Reached_Max {true};

        //Configs for prediction
        GenerateCandidateIsoform(RG,Candidate_Isfs,SubInsts);

        return;
    }

    void
    GenerateCandidateIsoform(format::ReadGroup& RG,TwoDimVec<bool>& Candidate_Isfs,std::vector<uint32_t>& SubInsts)
    {
        auto NumExons {RG.ExonBoundary.size()};
        auto MaxExonExplv {0.0},MaxJuncExonExplv {0.0};
        // HyperParameters
        double exon_min_frac         = EXON_MIN_FRAC , 
               junc_exp_frac         = JUNC_EXP_FRAC ,
               intron_retention_frac = INTRON_RETENTION_FRAC;
        std::vector<bool>   IsVisited(NumExons,false),
                            IntronRetention_Left(RG.ExonBoundary.size(),false),
                            IntronRetention_Right(RG.ExonBoundary.size(),false);
        std::vector<double> Explv,
                            JuncExplv,
                            LeftExplv(RG.ExonBoundary.size(),0),
                            RightExplv(RG.ExonBoundary.size(),0);
        std::vector<int32_t>  Indegree(RG.ExonBoundary.size(),0),
                              Outdegree(RG.ExonBoundary.size(),0);

                              
        //Calculate Expression levels for each exon.
        GetExpLv(RG,Explv,JuncExplv,LeftExplv,
                 RightExplv,Indegree,Outdegree,
                 IntronRetention_Left,IntronRetention_Right,
                 MaxExonExplv,MaxJuncExonExplv);

#ifdef DEBUG
        std::cout<<"Explv"<<std::endl;
        for(auto i:Explv)
            std::cout<<i<<" ";
        std::cout<<std::endl;

        std::cout<<"Junction level"<<std::endl;
        for(auto i:JuncExplv)
            std::cout<<i<<" ";
        std::cout<<std::endl;

        std::cout<<"Expression level"<<std::endl;
        for(auto i:Explv)
            std::cout<<i<<" ";
        std::cout<<std::endl;
#endif

        //IsoLasso::utils::Check_ReadGroup(RG);

        //Create tree structure
        std::vector<IsoLasso::Algorithm::ExonNode> ExonTree;
        for(auto Exon_idx=0;Exon_idx<NumExons;++Exon_idx)
            ExonTree.emplace_back(Exon_idx,RG.ReadSupportMatrix[Exon_idx],Explv,JuncExplv);


        //Choose different exons as beginning
        std::vector<bool>                  ExonsCoveredByIsoform(RG.ExonBoundary.size(),false);
        std::vector<std::vector<uint32_t>> Total_Paths;
        uint32_t        NumofFoundPaths       {0};

        /*
        for(auto i=0;i<ExonTree.size();i++)
        {
            std::cout<<"Exon "<<i<<std::endl;
            for(auto elem:ExonTree[i].getNeighbors())
                std::cout<<elem.first<<" ";
            std::cout<<'\n';
        }*/


        //Choose important exons and perform traversal
        while(true)
        {
            bool ReachedMax {false};
            std::vector<uint32_t> Selected_Exon;

            //Choose important exons as start using current parameters exon_min_frac,junc_exp_frac,
            for(auto exon_index=0;exon_index<RG.ExonBoundary.size();exon_index++)
            {
                bool Is_selected {((!(IntronRetention_Left[exon_index]||IntronRetention_Right[exon_index]) && Explv[exon_index]>=MaxExonExplv*exon_min_frac) // Is not intron retention
                                    ||((IntronRetention_Left[exon_index]||IntronRetention_Right[exon_index]) && Explv[exon_index]>=MaxExonExplv*intron_retention_frac) // Is intron retention
                                    ||(JuncExplv[exon_index] >= MaxJuncExonExplv*junc_exp_frac)
                                   ) && RG.ExonStats[exon_index][2] <= 0.7
                                 };
                
                if(Is_selected)
                    Selected_Exon.push_back(exon_index); // Should start from here
            }

            std::cout<<"SelectedExons"<<std::endl;
            IsoLasso::utils::print1Dvector(Selected_Exon);
            std::cout<<"ExonsCoveredByIsoform"<<std::endl;
            IsoLasso::utils::print1Dvector(ExonsCoveredByIsoform);

            bool FoundPath = false;
            for(auto Start:Selected_Exon) 
            {
                
                if(ExonsCoveredByIsoform[Start])
                    continue;
                else
                    FoundPath = true;
                NumofFoundPaths += TreeTraversal(Total_Paths,ExonTree,RG.ReadSupportMatrix,Explv,JuncExplv,Start,exon_min_frac,junc_exp_frac,intron_retention_frac,ExonsCoveredByIsoform);
                if(NumofFoundPaths>=MAX_ISOFORM_NUM)// Too many candidates
                {
                    Total_Paths.clear();
                    Total_Paths.shrink_to_fit();
                    break;
                }
                //IsoLasso::utils::print2Dvector(Total_Paths);
            }
            //Restrict Path 
            if(NumofFoundPaths>=MAX_ISOFORM_NUM)
            {
                NumofFoundPaths        = 0;
                exon_min_frac         *= THERSHOLD_GROWTH;
                intron_retention_frac *= THERSHOLD_GROWTH;

                if(intron_retention_frac<exon_min_frac)
                    intron_retention_frac = exon_min_frac;
                if(intron_retention_frac>1)
                    intron_retention_frac = 0.99;
                
                junc_exp_frac         *= THERSHOLD_GROWTH;
                if(junc_exp_frac>1)
                    junc_exp_frac      = 0.99;
            }
            if(!FoundPath)
                break;
         }
        
        std::cout<<"Found Paths:"<<std::endl;
        IsoLasso::utils::print2Dvector(Total_Paths);

        /*
        if(Total_Paths.size()!=0)
            FilterImpossibleCandidates(Candidate_Isfs,RG);
        */
        return;
   
    }

    uint32_t
    TreeTraversal(std::vector<std::vector<uint32_t>>& TotalPaths,
                  std::vector<IsoLasso::Algorithm::ExonNode>& Tree,
                  const std::vector<std::vector<uint32_t>>& ReadSupportMatrix,
                  const std::vector<double>& Explv,
                  const std::vector<double>& JuncExplv,
                  const uint32_t StartExon,
                  const double exon_min_frac,
                  const double junc_exp_frac ,
                  const double intron_retention_frac,
                  std::vector<bool>& ExonsCoveredByIsoform)
    {
        std::stack<std::pair<uint32_t,uint32_t>> ExonStack; // pair of {Exonindex,Depth}
        std::vector<uint32_t> CurrentPath;
        std::vector<bool> IsVisited(Tree.size(),false);
        uint32_t NumofPaths {0};

        ExonStack.emplace(StartExon,0);

        while(!ExonStack.empty())
        {
            //Get Node from top
            auto [CurrentExon,CurrentDepth] = ExonStack.top();
            ExonStack.pop();
            ExonsCoveredByIsoform[CurrentExon] = true;

            if(IsVisited[CurrentExon])
                continue;


            // Add new exon 
            CurrentPath.push_back(CurrentExon);
            IsVisited[CurrentExon] = true;            

            //Add neighbors to stack

            auto neighbors = Tree[CurrentExon].getNeighbors();

            auto ReachedEnd {true};

            if(!neighbors.empty())
            {
                auto MaxNeighborExplv  = std::ranges::max(neighbors| std::ranges::views::transform([](auto i){ return i.first; })
                                                                | std::ranges::views::transform([&Explv](auto i)->const auto&{ return Explv[i]; }));
                auto MaxNeighborJuncLv = neighbors.back().second;

                for(auto neighbor:neighbors)
                {            
                    if((Explv[neighbor.first]>=MaxNeighborExplv*exon_min_frac)&&
                    (ReadSupportMatrix[CurrentExon][neighbor.first]>=MaxNeighborJuncLv*junc_exp_frac)&&
                    !IsVisited[neighbor.first])
                    {
                        if(neighbor.first==(CurrentExon+1))//intron retention test
                        {
                            if(Explv[neighbor.first] < MaxNeighborExplv*intron_retention_frac)
                                continue;
                        }
                        ExonStack.emplace(neighbor.first,CurrentDepth+1);
                        ReachedEnd = false;
                    }
                    else    
                        continue;
                }
            }
            if(ReachedEnd) // end of path
            {
                TotalPaths.push_back(CurrentPath);
                for(auto Exon:CurrentPath)
                    ExonsCoveredByIsoform[Exon] = true;
                NumofPaths++;
                if(NumofPaths>=MAXPATH_PERNODE)
                {
                    std::cout<<"Warning: reached maximum path "<<MAXPATH_PERNODE <<" per node!"<<std::endl;
                    break;
                }
                else//Pop from stack
                {
                    if(ExonStack.empty())
                        break;

                    auto Pop_Count {CurrentDepth - ExonStack.top().second + 1};

                    while(Pop_Count--)
                    {
                        IsVisited[CurrentPath.back()] = false;
                        CurrentPath.pop_back();
                    }

                }
            }
        }
        return NumofPaths;
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
             std::vector<int32_t>& Indegree,
             std::vector<int32_t>& Outdegree,
             std::vector<bool>& IntronRetention_Left,
             std::vector<bool>& IntronRetention_Right,
             double& MaxExonExplv,
             double& MaxJuncExonExplv)
    {

        for(auto row=0;row<RG.ExonBoundary.size();++row)
        {
            for(auto col=row+1;col<RG.ExonBoundary.size();++col)
            {
                //From row exon to col exon
                if(RG.ReadSupportMatrix[row][col]>=MIN_JUNC_READ)
                {
                    Outdegree[row]++;
                    Indegree[col]++;
                }

                RightExplv[row] += RG.ReadSupportMatrix[row][col];
                LeftExplv[col]  += RG.ReadSupportMatrix[row][col];
            }
        }

        //Calculate each exon's expression level
        utils::GetExonExplv(RG,Explv);
        MaxExonExplv = std::ranges::max(Explv);

        //Junction
        for(auto exon_index=0;exon_index<RG.ExonBoundary.size();exon_index++)
        {
            LeftExplv[exon_index]  /= (RG.ReadLen*EXON_READLEN_FRAC);
            RightExplv[exon_index] /= (RG.ReadLen*EXON_READLEN_FRAC);
            JuncExplv.emplace_back((LeftExplv[exon_index]+RightExplv[exon_index])/2);
        }
        MaxJuncExonExplv = std::ranges::max(JuncExplv);

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
                       uint32_t& IsfStart, // Starting exon index
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
        Stack[0]      = IsfStart;
        Stack_Prev[0] = -1;

        while(Stackpt>=0) //Stack is not empty
        {
            uint32_t Current_exon {Stack[Stackpt]};
            bool     ReachedEnd   {true};
            Path[Current_exon]  =   true; //Visited

            if(!IsVisited[Stackpt];int32_t Prev_Stackpt=Stackpt)
            {
                IsVisited[Stackpt] = true;
                std::vector<uint32_t>& OutDegree = ReadSupportMatrix[Current_exon]; //OutDegree for current exon

                //Sort candidate exon by OutDegree value
                std::vector<uint32_t> ConnectedExons(OutDegree.size());
                for(auto exon_idx=0;exon_idx<OutDegree.size();++exon_idx)
                {
                    if(OutDegree[exon_idx]>=MIN_JUNC_READ)
                        ConnectedExons.emplace_back(exon_idx);
                }
                std::stable_sort(ConnectedExons.begin(), ConnectedExons.end(),
                                [OutDegree](auto i1, auto i2) {return OutDegree[i1] < OutDegree[i2];});

                if(!ConnectedExons.empty()) //Exist candidate exons
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
                            if(Explv[candidate_exon] >= MaxExpLv*EXON_MIN_FRAC)
                                ValidConnectedExons[index] = true;
                            else
                                ValidConnectedExons[index] = false;
                        }
                    }
                    //IntronRetention
                    auto NextExon {std::find(ConnectedExons.begin(),ConnectedExons.end(),Current_exon+1)};
                    if(NextExon != ConnectedExons.end()) //Intron retention exist
                    {
                        if(Explv[*NextExon]<MaxExpLv*intronretentionfrac)
                            ValidConnectedExons[std::distance(ConnectedExons.begin(),NextExon)] = false;
                    }

                    //Filter exons with ValidConnectedExons = false
                    std::erase_if(ConnectedExons,[&ValidConnectedExons](auto index){return !ValidConnectedExons[index];});
                    
                    ReachedEnd = ConnectedExons.empty();

                    for(auto index:ConnectedExons)//Put all valid connected exons into stack
                    {
                        Stack.emplace_back(index);
                        Stack_Prev.emplace_back(Prev_Stackpt);
                        IsVisited[index] = false;
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
                    if(NumIsofs>=MAX_ISOFORM_NUM)
                        break;
                }
            }
            else
                Stackpt--;
        }
        return;
    }

    /*
     * 
     *
     */

    void
    FilterImpossibleCandidates(TwoDimVec<bool>& Candidaite_Isfs,format::ReadGroup& RG)
    {
        auto NumExons {RG.ExonBoundary.size()},NumCandidates {Candidaite_Isfs.size()};
        std::vector<bool> IsValid(NumCandidates,false);

        std::vector<std::vector<uint32_t>> SGJuncType;
        std::vector<uint32_t>              SGJuncCount;

        for(auto Typeindex=0;Typeindex<RG.SGTypes.size();Typeindex++)
        {
            if(std::accumulate(RG.SGTypes[Typeindex].begin(),
                               RG.SGTypes[Typeindex].end(),
                               0,
                               [](auto i1,auto i2){return i1+(i2==true);})
                               >1)//With more than one exon
            {
                SGJuncType.emplace_back(RG.SGTypes[Typeindex]);
                SGJuncCount.emplace_back(RG.TypeCount[Typeindex]);
            }
        }

        for(auto Candidateindex=0;Candidateindex<NumCandidates;Candidateindex++)
        {
            if(std::accumulate(Candidaite_Isfs[Candidateindex].begin(),
                               Candidaite_Isfs[Candidateindex].end(),
                               0,
                               [](auto i1,auto i2){return i1+(i2==true);})
                               ==1)//Single exon
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