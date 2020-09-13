#include <SAM_module/ReadGroup.hpp>
#include <utils/Commontype.hpp>
#include <SAM_module/RangeSet.hpp>
#include <utils/Auxiliario.hpp>

std::ofstream IsoLasso::utils::RG_STATS_FS;


namespace IsoLasso::format
{
    void
    ReadGroup::AddRecord(const IsoLasso::format::Sam_record& record,range_type& current_range)
    {
        AddWithoutPair(record);
        if(record.isPairedEnd)
        {
            if(record.RNext==record.RName || record.RNext=="=")
            {
                if(QNameQueryTable[record.QName].find(record.Pos)!=QNameQueryTable[record.QName].end())
                {
                    PairendTable[QNameQueryTable[record.QName][record.Pos]]=ReadStart.size()-1;
                    PairendTable.back() = QNameQueryTable[record.QName][record.Pos];
                    QNameQueryTable[record.QName].erase(QNameQueryTable[record.QName][record.Pos]);                    
                }
                else
                    QNameQueryTable[record.QName][record.PNext] = ReadStart.size()-1;
            }
        }
        return;
    }
    void
    ReadGroup::RemoveLongSpanReads(const uint32_t& MAX_EXON_SPAN,const uint32_t& MAX_JUNCTION_COVERAGE)
    {
        std::map<u_int32_t,u_int32_t> Coverage;
        for(auto read:ReadStart) 
            Coverage[read[0]]++;
        for(auto read_idx = 0;read_idx<ReadStart.size();read_idx++)
        {
            if(ValidRead[read_idx])
            {
                for(auto exon_idx=1;exon_idx<ReadStart[read_idx].size();exon_idx++)
                {
                    u_int32_t cur_start{ReadStart[read_idx][exon_idx]}, prev_end {ReadEnd[read_idx][exon_idx-1]};
                    u_int32_t span {cur_start-prev_end};
                    if(span<MAX_EXON_SPAN)
                        continue;
                    else if(utils::ExceedThreshold(Coverage,prev_end,cur_start,MAX_JUNCTION_COVERAGE))
                    {
                        ValidRead[read_idx] = false;
                        if(PairendTable[read_idx]!=-1)
                            ValidRead[PairendTable[read_idx]] = false;
                        break;
                    }
                }
            }
            if(PairendTable[read_idx]!=-1&&read_idx<PairendTable[read_idx])
            {
                u_int32_t left_end_back {ReadEnd[read_idx].back()},right_end_front {ReadStart[PairendTable[read_idx]].front()};
                u_int32_t pair_span {right_end_front-left_end_back};
                if(pair_span<MAX_EXON_SPAN)
                    continue;
                else if(utils::ExceedThreshold(Coverage,left_end_back,right_end_front,MAX_JUNCTION_COVERAGE))
                {
                    ValidRead[read_idx] = false;
                    ValidRead[PairendTable[read_idx]] = false;
                }
            }
        }
        return;
    }
    void
    ReadGroup::SplitbyRangeSet(std::vector<ReadGroup>& SubRGs,const uint32_t& MIN_GAP_SPAN)
    {
        std::vector<IsoLasso::format::RangeSet> RangeSets;
        for(auto read_index=0;read_index<ReadStart.size();read_index++)
        {
            if(ValidRead[read_index])
            {
                bool NewRange {true};
                for(auto range_index=0;range_index<RangeSets.size();range_index++)
                {
                    if(PairendTable[read_index]==-1)
                    {
                        u_int32_t current_minDist {RangeSets[range_index].MinDistance(ReadStart[read_index],ReadEnd[read_index])};
                        if(current_minDist<= MIN_GAP_SPAN)
                        {
                            RangeSets[range_index].Add(ReadStart[read_index],ReadEnd[read_index]);
                            SubRGs[range_index].AddWithoutPair(ReadStart[read_index],ReadEnd[read_index],Direction[read_index]);
                            NewRange = false;
                            break;
                        }
                    }
                    else
                    {
                        int32_t Paired_index {PairendTable[read_index]};
                        if(read_index<Paired_index)
                        {
                            if(RangeSets[range_index].MinDistance(ReadStart[read_index],ReadEnd[read_index])<=MIN_GAP_SPAN ||
                               RangeSets[range_index].MinDistance(ReadStart[Paired_index],ReadEnd[Paired_index]))
                            {
                                RangeSets[range_index].Add(ReadStart[read_index],ReadEnd[read_index]);
                                RangeSets[range_index].Add(ReadStart[Paired_index],ReadEnd[Paired_index]);
                                SubRGs[range_index].AddPair(ReadStart[read_index],ReadEnd[read_index],Direction[read_index],
                                                            ReadStart[Paired_index],ReadEnd[Paired_index],Direction[Paired_index]);
                                NewRange = false;
                                break;
                            }
                        }
                    }
                }
                if(NewRange)
                {
                    if(PairendTable[read_index]==-1 || read_index<PairendTable[read_index])
                    {
                        ReadGroup NewRG;
                        NewRG.ChrName = ChrName;
                        RangeSet NewRangeSet;
                        NewRG.RG_index = RG_index;

                        if(PairendTable[read_index]!=-1)
                        {
                            int32_t Paired_index {PairendTable[read_index]};
                            NewRangeSet.Add(ReadStart[read_index],ReadEnd[read_index]);
                            NewRangeSet.Add(ReadStart[Paired_index],ReadEnd[Paired_index]);
                            NewRG.AddPair(ReadStart[read_index],ReadEnd[read_index],Direction[read_index],
                                          ReadStart[Paired_index],ReadEnd[Paired_index],Direction[Paired_index]);
                        }
                        else
                        {
                            NewRangeSet.Add(ReadStart[read_index],ReadEnd[read_index]);
                            NewRG.AddWithoutPair(ReadStart[read_index],ReadEnd[read_index],Direction[read_index]);
                        }
                        SubRGs.emplace_back(NewRG);
                        RangeSets.emplace_back(NewRangeSet);
                    }
                }
            }
        }
        std::vector<bool> ValidRanges(RangeSets.size(),true);
        bool Merged {true};
        while(Merged)
        {
            Merged = false;
            for (auto outer=0;outer<RangeSets.size();outer++)
            {
                if(ValidRanges[outer])
                {
                    for(auto inner=outer+1;inner<RangeSets.size();inner++)
                    {
                        if(ValidRanges[inner])
                        {
                            if(RangeSets[outer].MinDistance(RangeSets[inner])<=MIN_GAP_SPAN)
                            {
                                auto [InnerStart,InnerEnd] {RangeSets[inner].toRange()};
                                RangeSets[outer].Add(InnerStart,InnerEnd);
                                for(auto read_index=0;read_index<SubRGs[inner].ReadStart.size();read_index++)
                                {
                                    int32_t paired_index {SubRGs[inner].PairendTable[read_index]};
                                    if(paired_index==-1)
                                        SubRGs[outer].AddWithoutPair(SubRGs[inner].ReadStart[read_index],SubRGs[inner].ReadEnd[read_index],Direction[read_index]);
                                    else
                                    {
                                        if(read_index<paired_index)
                                            SubRGs[outer].AddPair(SubRGs[inner].ReadStart[read_index],SubRGs[inner].ReadEnd[read_index],Direction[read_index],
                                                                  SubRGs[inner].ReadStart[paired_index],SubRGs[inner].ReadEnd[paired_index],Direction[paired_index]);
                                    }
                                }
                                ValidRanges[inner] = false;
                                Merged = true;
                            }
                        }
                    }
                }
            }
        }
        std::vector<ReadGroup>::iterator RG_iter {SubRGs.begin()};
        uint32_t Sub_index {0};

        for(auto RG_index =0;RG_index<ValidRanges.size();RG_index++)
        {
            if(ValidRanges[RG_index] && RG_iter->ValidRead.size()>MIN_RG_SIZE)
            {
                //Initial orientation
                RG_iter->SetOrientation();   
                RG_iter->CurrentRange=RG_iter->getRange();
                RG_iter->SubRG_index = Sub_index;
                Sub_index++;
                RG_iter++;
            }
            else
                RG_iter = SubRGs.erase(RG_iter);
        }
        return;
    }
    void
    ReadGroup::SplitbyDirection(std::vector<ReadGroup>& Split_SubRGs,std::vector<range_type>& Split_Ranges)
    {
        //TODO
    }
    void 
    ReadGroup::CalculateBound(const uint32_t MIN_JUNC_COV,const uint32_t MIN_GAP_SPAN)
    {
        std::map<uint32_t,uint32_t> SegBoundary;
        for(auto read_index=0;read_index<ReadStart.size();read_index++)
        {
            if(ReadStart[read_index].size()>1)// Only junction reads can be used to calculate boundary
            {
                for(auto exon_index=1;exon_index<ReadStart[read_index].size()-1;exon_index++)
                {
                    SegBoundary[ReadStart[read_index][exon_index]]++;
                    SegBoundary[ReadEnd[read_index][exon_index]+1]++;
                }
                SegBoundary[ReadStart[read_index].back()]++;
                SegBoundary[ReadEnd[read_index].front()+1]++;                
            }
        }
        std::map<uint32_t,uint32_t> Boundary;
        Boundary[CurrentRange.first] = 0;
        Boundary[CurrentRange.second+1] = 0;
        for(auto map_iter=SegBoundary.begin();map_iter!=SegBoundary.end();map_iter++)
        {
            if(map_iter->second>=MIN_JUNC_COV)
                Boundary[map_iter->first] = 0;
        }
        SegBoundary.clear();
        std::map<uint32_t,int32_t> coverage;
        GetCoverage(coverage);
        std::vector<range_type> Cutpoint;
        GetCvgCutPoint(coverage,Cutpoint,MIN_EXON_COV,MIN_GAP_SPAN);
        for(auto cvg_range:Cutpoint)
        {
            if(Boundary.count(cvg_range.first)==0)
                Boundary[cvg_range.first]++;
            if(Boundary.count(cvg_range.second+1)==0)
                Boundary[cvg_range.second+1]++;
        }
        GetCvgStats(coverage,Boundary,CvgStats);
        //Save boundary as rangetype
        std::map<u_int32_t,u_int32_t>::const_iterator Bound_iter {Boundary.begin()},
                                                      next_Bound_iter {std::next(Bound_iter)};
        while(next_Bound_iter!=Boundary.end())
        {
            ExonBoundary.emplace_back(range_type(Bound_iter->first,next_Bound_iter->first-1));
            Bound_iter++;
            next_Bound_iter++;
        }
        ValidExons.assign(ExonBoundary.size(),true);
        return;
    }
    void
    ReadGroup::GetCoverage(std::map<uint32_t,int32_t>& coverage)
    {
        //Places to be calculated
        for(auto read_index=0;read_index<ReadStart.size();read_index++)
        {
            if(ValidRead[read_index])
            {
                for(auto exon_index=0;exon_index<ReadStart[read_index].size();exon_index++)
                {
                    coverage[ReadStart[read_index][exon_index]]+=1;
                    coverage[ReadEnd[read_index][exon_index]+1]-=1;
                }
            }
        }
        //Calculate coverage w/ ++-- algirithm
        int32_t cvg_count {0};
        for(auto cvg_iter=coverage.begin();cvg_iter!=coverage.end();)
        {
            if(cvg_iter->second!=0)
            {
                cvg_count +=cvg_iter->second;
                cvg_iter->second = cvg_count;
                cvg_iter++;
            }
            else //Condense
                coverage.erase(cvg_iter++);
        }
        return;
    }
    void 
    ReadGroup::GetCvgCutPoint(const std::map<uint32_t,int32_t>&coverage,std::vector<range_type>& Cutpoint,
                   const u_int32_t threshold,const u_int32_t MIN_GAP_SPAN)
    {
        std::vector<range_type> CvgRanges;
        bool NewRange {false};
        for(auto cvg_iter=coverage.begin();cvg_iter!=coverage.end();cvg_iter++)
        {
            if(NewRange)
            {
                if(cvg_iter->second<threshold)
                    NewRange = false;
                else //extend
                    CvgRanges.back().second = cvg_iter->first;
            }
            else
            {
                if(cvg_iter->second>=threshold)
                {
                    CvgRanges.emplace_back(cvg_iter->first,cvg_iter->first+1);
                    NewRange = true;
                }
            }
        }
        for(auto Cur_range:CvgRanges)
        {
            if(Cutpoint.size()>0 && Cur_range.first-Cutpoint.back().second <=MIN_GAP_SPAN)//Merge
                Cutpoint.back().second = Cur_range.second;
            else//NewRange
                Cutpoint.emplace_back(Cur_range);
        }
        return;
    }
    void 
    ReadGroup::GetCvgStats( std::map<uint32_t,int32_t>&coverage,const std::map<uint32_t,uint32_t>& Boundary,
                            std::map<uint32_t,std::vector<double>>& CvgStats)
    {
        //Merging Boundary into coverage
        for(auto bound_iter=Boundary.begin();bound_iter!=Boundary.end();bound_iter++)
        {
            if(coverage.count(bound_iter->first))
                continue;
            else
            {
                auto insert {coverage.insert(std::pair<uint32_t,int32_t>(bound_iter->first,0))};
                if(insert.second==false)//Exist
                    continue;
                else if(insert.first->first!=(coverage.begin())->first) // Add to coverage
                    insert.first->second = std::prev(insert.first)->second;
            }
        }
        auto Bound_iter  {Boundary.begin()};
        auto next_Bound_iter {std::next(Bound_iter)};
        //Calculate status for every exon
        while(next_Bound_iter!=Boundary.end())
        {
            auto Cvg_iter  {coverage.find(Bound_iter->first)};
            auto next_cvg_iter {coverage.find(next_Bound_iter->first)};
            auto Inneriter {Cvg_iter},Inneriter2 {std::next(Cvg_iter)};
            double  StartCvg    {double(Cvg_iter->second)}, 
                    Zerofrac    {0}, 
                    MaxCvg      {-1},
                    MeanCvg     {0};
            //Cvg stats within an exon
            while(Inneriter!=next_cvg_iter)
            {
                if(Inneriter->second>MaxCvg)
                    MaxCvg = Inneriter->second; 
                if(Inneriter->second==0)//length of zero coverage
                    Zerofrac+=(Inneriter2->first-Inneriter->first);
                MeanCvg+=(Inneriter2->first-Inneriter->first)*Inneriter->second;
                Inneriter++;
                Inneriter2++;
                if(Inneriter2==coverage.end())
                    break;
            }
            CvgStats[Bound_iter->first]=std::vector<double>{MaxCvg,
                                                            StartCvg,
                                                            Zerofrac/(next_cvg_iter->first-Cvg_iter->first),
                                                            MeanCvg/(next_cvg_iter->first-Cvg_iter->first)};
            Bound_iter++;
            next_Bound_iter++;
        }
        return;
    }
    void
    ReadGroup::CalculateType()
    {
        ExonCoverage.assign(ExonBoundary.size(),0);
        for(auto read_index=0;read_index<ReadStart.size();read_index++)
        {
            std::vector<uint32_t> ReadType {GetType(ReadStart[read_index],ReadEnd[read_index],ExonCoverage,MIN_OVERLAP)};
            auto TypeIter {std::find(SGTypes.begin(),SGTypes.end(),ReadType)};
            if(TypeIter==SGTypes.end())
            {
                Read2Type.emplace_back(SGTypes.size());
                SGTypes.emplace_back(ReadType);                
                TypeCount.emplace_back(1);
                TypeDirection.emplace_back(Direction[read_index]);
                ValidType.emplace_back(ReadType.size()>0);
            }
            else
            {
                auto TypeIndex {std::distance(SGTypes.begin(),TypeIter)};
                TypeCount[TypeIndex]++;
                TypeDirection[TypeIndex]+=Direction[read_index];
                Read2Type.emplace_back(TypeIndex);
            }
        }
        return;
    }
    std::vector<uint32_t>
    ReadGroup::GetType(const std::vector<uint32_t>& SegStart,const std::vector<uint32_t>& SegEnd,
                        std::vector<uint32_t>& ExonCoverage,const uint32_t& MIN_OVERLAP)
    {
        std::vector<uint32_t> ReadType;
        for(auto seg_index=0;seg_index<SegStart.size();seg_index++)
        {
            for(auto exon_index=0;exon_index<ExonBoundary.size();exon_index++)
            {
                uint32_t OverlapDist {utils::GetOverLapping(SegStart[seg_index],SegEnd[seg_index],
                                                               ExonBoundary[exon_index].first,ExonBoundary[exon_index].second)};
                if(OverlapDist>0)//Is overlapping
                {
                    if(OverlapDist<MIN_OVERLAP)//Not Large enough
                    {
                        //3' site
                        if( seg_index==0 &&
                            SegEnd[seg_index]>ExonBoundary[exon_index].second &&
                            ExonBoundary[exon_index].second>SegStart[seg_index] &&
                            SegStart[seg_index]>ExonBoundary[exon_index].first)
                            continue;
                        if( seg_index==SegStart.size()-1 &&
                            SegEnd[seg_index]<ExonBoundary[exon_index].second &&
                            ExonBoundary[exon_index].first<SegEnd[seg_index] &&
                            SegStart[seg_index]<ExonBoundary[exon_index].first)
                            continue;
                    }
                }
                ExonCoverage[exon_index]++;
                ReadType.emplace_back(exon_index);
            }
        }
        std::sort( ReadType.begin(), ReadType.end() );
        ReadType.erase( std::unique( ReadType.begin(), ReadType.end() ), ReadType.end() );        
        return ReadType;
    }
    void
    ReadGroup::RemoveWeakExons(const double MIN_CVG_FRAC)
    {
        std::vector<bool> Include_InJunction(ExonBoundary.size(),0);

        for(auto Type:SGTypes)
        {
            if(Type.size()>1)
            {
                if(Type[0]!=Type[1]-1)
                    Include_InJunction[Type[0]] = true;
                if(Type.back()!=*std::prev(Type.end())+1)
                    Include_InJunction[Type.back()] = true;
                for(auto exon_index=1;exon_index<Type.size()-1;exon_index++)
                    Include_InJunction[Type[exon_index]] = true;
            }
        }

        for(auto exon_index=0;exon_index<ExonBoundary.size();exon_index++)
        {
            if(ExonCoverage[exon_index]==0)
            {
                ValidExons[exon_index] = false;
                continue;
            }

            double Neighbor_Cvg {0.0};

            if(exon_index==0)
                Neighbor_Cvg = CvgStats[ExonBoundary[1].first][3];
            else if(exon_index==ExonBoundary.size()-1)
                Neighbor_Cvg = CvgStats[ExonBoundary[ExonBoundary.size()-2].first][3];
            else
                Neighbor_Cvg = CvgStats[ExonBoundary[exon_index+1].first][3]>CvgStats[ExonBoundary[exon_index-1].first][3]?
                               CvgStats[ExonBoundary[exon_index+1].first][3]:
                               CvgStats[ExonBoundary[exon_index-1].first][3];

            if(CvgStats[ExonBoundary[exon_index].first][3]<=Neighbor_Cvg*MIN_CVG_FRAC 
               && !Include_InJunction[exon_index])
                ValidExons[exon_index] = false; 
        }
        return;
    }
    void
    ReadGroup::CalculateValidExons()
    {
        bool HasChanged {true};

        while(HasChanged)
        {
            for(auto Type_index=0;Type_index<SGTypes.size();Type_index++)
            {
                for(auto exon_index=0;exon_index<SGTypes[Type_index].size();exon_index++)
                {
                    if(!ValidExons[SGTypes[Type_index][exon_index]])
                    {
                        ValidType[Type_index] = false;
                        break;
                    }
                }
            }

            ExonCoverage.assign(ExonCoverage.size(),0);
            HasChanged = false;

            for(auto read_index=0;read_index<ReadStart.size();read_index++)
            {
                if(!ValidType[Read2Type[read_index]])
                {
                    ValidRead[read_index] = false;
                    continue;
                }
                for(auto exon_index=0;exon_index<SGTypes[Read2Type[read_index]].size();exon_index++)
                    ExonCoverage[SGTypes[Read2Type[read_index]][exon_index]]++;
            }

            for(auto exon_index=0;exon_index<ExonBoundary.size();exon_index++)
            {
                if(!ValidExons[exon_index])
                    continue;
                else if(ExonCoverage[exon_index]==0)
                {
                    ValidExons[exon_index] = false;
                    HasChanged = true;
                }
            }
        }

        std::map<uint32_t,int32_t> Coverage;
        GetCoverage(Coverage);

        CvgStats.clear();
        std::map<uint32_t,uint32_t> Boundary;

        for(auto bound:ExonBoundary)
        {
            Boundary[bound.first]++;
            Boundary[bound.second+1]++;
        }
        GetCvgStats(Coverage,Boundary,CvgStats);
        return;
    }

    void
    ReadGroup::WriteStatsToFile(std::ofstream& ofs)
    {
        ofs<<"ReadGroup "<<RG_index+1<<"-"<<SubRG_index+1<<std::endl;
        ofs<<"Boundary "<<ChrName<<"\t"<<CurrentRange.first<<" "<<CurrentRange.second<<" "<<Orientation<<std::endl;
        ofs<<"Number of exons:"<<ExonBoundary.size()<<std::endl;
        ofs<<"Total valid reads:"<<validSize()<<std::endl;
        for(auto exon_index=0;exon_index<ExonBoundary.size();exon_index++)
        {
            if(ValidExons[exon_index])
            {
                ofs<<ExonBoundary[exon_index].first<<" "<<ExonBoundary[exon_index].second<<" "<<(ExonBoundary[exon_index].second-ExonBoundary[exon_index].first+1)<<"　"<<ExonCoverage[exon_index]<<" ";
                auto Cvg_iter {CvgStats.find(ExonBoundary[exon_index].first)};
                for(auto stats:Cvg_iter->second)
                    ofs<<stats<<",";
                ofs<<std::endl;
            }
        }

        ofs<<"SGTypes:"<<std::accumulate(ValidType.begin(),ValidType.end(),0,[](uint32_t a,bool b){return a+(b==true);})<<std::endl;
        for(auto TypeIndex=0;TypeIndex<SGTypes.size();TypeIndex++)
        {
            if(!ValidType[TypeIndex])
                continue;

            auto Exon_iter {SGTypes[TypeIndex].begin()};

            for(auto exon_index=0;exon_index<ExonBoundary.size();exon_index++)
            {
                if(ValidExons[exon_index])
                {
                    if(*Exon_iter==exon_index)
                    {
                        ofs<<1<<" ";
                        Exon_iter++;
                    }
                    else
                        ofs<<0<<" ";
                }
            }
            ofs<<TypeCount[TypeIndex]<<" ";
            if(TypeDirection[TypeIndex]>0)
                ofs<<"+"<<std::endl;
            else if(TypeDirection[TypeIndex]<0)
                ofs<<"-"<<std::endl;
            else
                ofs<<"."<<std::endl;
        }

        uint32_t nPETypes {0};
        std::map<range_type,std::vector<uint32_t>> PETypes;

        for(auto read_index=0;read_index<ReadStart.size();read_index++)
        {
            if(PairendTable[read_index]!=-1 && read_index<PairendTable[read_index])
            {
                uint32_t LeftType  {Read2Type[read_index]};
                uint32_t RightType {Read2Type[PairendTable[read_index]]};

                if(!(ValidType[LeftType]&&ValidType[RightType]))
                    continue;

                nPETypes++;
                PETypes[range_type(LeftType,RightType)].emplace_back(ReadStart[PairendTable[read_index]].front()-ReadEnd[read_index].back());
            }
        }

        uint32_t TypeCount {0};
        std::vector<uint32_t> ValidSGIndex(ValidType.size(),0);
        for(auto Type_index=0;Type_index<ValidType.size();Type_index++)
        {
            if(ValidType[Type_index])
                ValidSGIndex[Type_index] = TypeCount++;
            else
                ValidSGIndex[Type_index] = -1;
        }


        ofs<<"PETypes:"<<nPETypes<<std::endl;
        for(auto PE_Iter=PETypes.begin();PE_Iter!=PETypes.end();PE_Iter++)
        {
            ofs<<ValidSGIndex[PE_Iter->first.first]<<" "<<
                 ValidSGIndex[PE_Iter->first.second]<<" "<<
                 PE_Iter->second.size()<<" ";
            ofs<<std::accumulate(PE_Iter->second.begin(),PE_Iter->second.end(),0.0)/PE_Iter->second.size()<<std::endl;
        }

        ofs<<"Coverage for SGTypes:"<<std::endl;

        




        return;
    }   

}//end of IsoLasso::format
namespace IsoLasso::utils
{
    void
    ProcessReadGroup(IsoLasso::format::ReadGroup& RG,const range_type& current_range)
    {
        std::vector<format::ReadGroup> SubRGs;
        RG.RemoveLongSpanReads(MAX_EXON_SPAN,MAX_JUNCTION_COVERAGE);
        RG.SplitbyRangeSet(SubRGs,MIN_GAP_SPAN);
        for(auto SubRG:SubRGs)
        {
            //Enumerate all boundaries
            SubRG.CalculateBound(MIN_JUNC_COV,MIN_GAP_SPAN);
            //Calculate SGType (Combinations of high quality exon combinations)
            SubRG.CalculateType();
            //Remove Exons with too small coverage and not included in any junction type
            SubRG.RemoveWeakExons(MIN_CVG_FRAC);
            //Remove Types with Invalid exons and those reads with invalid types and recalculate coverage statistics
            SubRG.CalculateValidExons();

            if(SubRG.validSize()<MIN_RG_SIZE)
                continue;

            SubRG.WriteStatsToFile(RG_STATS_FS);
        }
        return ;
    }


    

}//end of IsoLasso::utils