#include <SAM_module/ReadGroup.hpp>
#include <utils/Commontype.hpp>
#include <SAM_module/RangeSet.hpp>
#include <utils/Auxiliario.hpp>
#include <SAM_module/SAMrecord.hpp>
#include <EM_module/PredExpLevel.hpp>
#include <EM_module/EMalgorithm.hpp>
#include <utils/ThreadPool.hpp>
#include <algorithm>

namespace IsoLasso::format
{
    void
    ReadGroup::AddRecord(const IsoLasso::format::Sam_record& record)
    {
        AddWithoutPair(record);
        
        if(record.isPairedEnd)
        {        
            if(QNameQueryTable[record.QName].find(record.Pos)!=QNameQueryTable[record.QName].end())//Other hand not exist
            {
                PairendTable[QNameQueryTable[record.QName][record.Pos]] = ReadStart.size()-1; //Current Read
                PairendTable.back() = QNameQueryTable[record.QName][record.Pos];
                QNameQueryTable[record.QName].erase(QNameQueryTable[record.QName][record.Pos]);                    
            }
            else
                QNameQueryTable[record.QName][record.PNext] = ReadStart.size()-1; //Current Read 
        }
        ReadLen_Count[utils::GetEfficientLen(record)]++;
        return;
    }
    void
    ReadGroup::RemoveLongSpanReads(const uint32_t& MAX_EXON_SPAN,const uint32_t& MAX_JUNCTION_COVERAGE)
    {
        std::map<uint32_t,uint32_t> Coverage;
        for(const auto& read:ReadStart) 
            Coverage[read[0]]++;
        for(auto read_idx = 0;read_idx<ReadStart.size();read_idx++)
        {
            if(ValidRead[read_idx])
            {
                for(auto seg_idx=1;seg_idx<ReadStart[read_idx].size();seg_idx++)//Segments
                {
                    uint32_t cur_start{ReadStart[read_idx][seg_idx]}, prev_end {ReadEnd[read_idx][seg_idx-1]}; 
                    uint32_t span {cur_start-prev_end};// Spans for two segments
                    if(span<MAX_EXON_SPAN)
                        continue;
                    else if(utils::ExceedThreshold(Coverage,prev_end,cur_start,MAX_JUNCTION_COVERAGE))
                    {
                        ValidRead[read_idx] = false;
                        if(PairendTable[read_idx]!=-1) // Paired-End
                            ValidRead[PairendTable[read_idx]] = false;
                        break;
                    }
                }
            }
            if(PairendTable[read_idx]!=-1&&read_idx<PairendTable[read_idx]) // Paired-end's left side
            {
                uint32_t left_end_back {ReadEnd[read_idx].back()},right_end_front {ReadStart[PairendTable[read_idx]].front()};
                uint32_t pair_span {right_end_front-left_end_back};
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
        for(auto read_index=0;read_index<ReadStart.size();read_index++) //For each read
        {
            if(ValidRead[read_index])
            {
                bool NewRange {true};
                for(auto range_index=0;range_index<RangeSets.size();range_index++)//For every existing range
                {
                    if(PairendTable[read_index]==-1) // Single-End
                    {
                        uint32_t current_minDist {RangeSets[range_index].MinDistance(ReadStart[read_index],ReadEnd[read_index])};
                                                
                        if(current_minDist<= MIN_GAP_SPAN)//Redistribute reads
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
                               RangeSets[range_index].MinDistance(ReadStart[Paired_index],ReadEnd[Paired_index])<=MIN_GAP_SPAN)
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
                        SubRGs.push_back(NewRG);
                        RangeSets.push_back(NewRangeSet);
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
                                auto [InnerStart,InnerEnd] {std::move(RangeSets[inner].toRange())};
                                RangeSets[outer].Add(InnerStart,InnerEnd);
                                for(auto read_index=0;read_index<SubRGs[inner].ReadStart.size();read_index++)
                                {
                                    int32_t paired_index {SubRGs[inner].PairendTable[read_index]};
                                    if(paired_index==-1)
                                        SubRGs[outer].AddWithoutPair(SubRGs[inner].ReadStart[read_index],
                                                                     SubRGs[inner].ReadEnd[read_index],
                                                                     Direction[read_index]);
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

        for(auto RG_Index =0;RG_Index<ValidRanges.size();RG_Index++)
        {
            if(ValidRanges[RG_Index] && RG_iter->ValidRead.size()>MIN_RG_SIZE)
            {
                //Initial orientation
                RG_iter->SetOrientation();   
                RG_iter->CurrentRange = std::move(RG_iter->getRange());
                RG_iter->SubRG_index  = Sub_index;
                RG_iter->ReadLen = std::max_element(ReadLen_Count.begin(), ReadLen_Count.end(),
                                                    [](const std::pair<uint32_t, uint32_t>& p1, const std::pair<uint32_t, uint32_t>& p2) {
                                                    return p1.second < p2.second; })->first;
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
        //Get possible boundaries from segments.
        std::map<uint32_t,uint32_t> SegBoundary;
        for(auto read_index=0;read_index<ReadStart.size();read_index++)//Every read
        {
            if(ReadStart[read_index].size()>1)// Only junction reads can be used to calculate boundary
            {
                for(auto exon_index=1;exon_index<ReadStart[read_index].size()-1;exon_index++)
                {
                    SegBoundary[ReadStart[read_index][exon_index]]++;
                    SegBoundary[ReadEnd[read_index][exon_index]+1]++;//ReadEnd[read_index][exon_index] is the last aligned base 
                }
                SegBoundary[ReadStart[read_index].back()]++;
                SegBoundary[ReadEnd[read_index].front()+1]++;                
            }
        }

        std::map<uint32_t,uint32_t> Boundary;
        Boundary[CurrentRange.first] = 0;
        Boundary[CurrentRange.second+1] = 0;

        //Set boundaries with high coverage
        for(auto map_iter=SegBoundary.begin();map_iter!=SegBoundary.end();map_iter++)
        {
            if(map_iter->second>=MIN_JUNC_COV)
                Boundary[map_iter->first] = 0; // Add Boundary
        }
        SegBoundary.clear();

        //Get coverage by ++-- algorithm
        std::map<uint32_t,int32_t> coverage;
        GetCoverage(coverage);

        //Get boundary by coverage cutpoint
        std::vector<range_type> Cutpoint;
        GetCvgCutPoint(coverage,Cutpoint,MIN_EXON_COV,MIN_GAP_SPAN);

        // Update Coverage Cutpoint to boundary
        // Do not insert boundaries if the maximum coverage is too low, or if the gap is too short
        // Check adjacent ranges, if they're too close, merge them
        for(const auto& cvg_range:Cutpoint)
        {
            if(Boundary.count(cvg_range.first)==0 && utils::ShortestDist(Boundary,cvg_range.first)>ReadLen)
                Boundary[cvg_range.first]++;
            if(Boundary.count(cvg_range.second+1)==0 && utils::ShortestDist(Boundary,cvg_range.second)>ReadLen)
                Boundary[cvg_range.second+1]++;
        }

        GetCvgStats(coverage,Boundary,CvgStats);
        //Save boundary as rangetype
        auto Bound_iter {Boundary.begin()},
             next_Bound_iter {std::next(Bound_iter)};
        while(next_Bound_iter!=Boundary.end())
        {
            ExonBoundary.emplace_back(Bound_iter->first,next_Bound_iter->first-1);
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
                for(auto seg_index=0;seg_index<ReadStart[read_index].size();seg_index++) //Total segments for current read
                {
                    coverage[ReadStart[read_index][seg_index]] += 1;
                    coverage[ReadEnd[read_index][seg_index]+1] -= 1;
                }
            }
        }

        //Calculate coverage w/ ++-- algirithm
        int32_t cvg_count {0};
        for(auto cvg_iter=coverage.begin();cvg_iter!=coverage.end();)
        {
            if(cvg_iter->second!=0)
            {
                cvg_count += cvg_iter->second;
                cvg_iter->second = cvg_count;
                cvg_iter++;
            }
            else //Condense
                coverage.erase(cvg_iter++);
        }
        return;
    }
    void 
    ReadGroup::GetCvgCutPoint(const std::map<uint32_t,int32_t>&coverage,
                                    std::vector<range_type>& Cutpoint,
                              const uint32_t MIN_EXON_COV,
                              const uint32_t MIN_GAP_SPAN)
    {
        std::vector<range_type> CvgRanges;
        bool NewRange {false};
        for(auto cvg_iter=coverage.begin();cvg_iter!=coverage.end();cvg_iter++)
        {
            if(NewRange)
            {
                if(cvg_iter->second<MIN_EXON_COV) 
                    NewRange = false;
                else //extend range
                    CvgRanges.back().second = cvg_iter->first;
            }
            else
            {
                if(cvg_iter->second>=MIN_EXON_COV) // Start a New Range
                {
                    CvgRanges.emplace_back(cvg_iter->first,cvg_iter->first+1);
                    NewRange = true;
                }
            }
        }

        // Merge 
        for(const auto& Cur_range:CvgRanges)
        {
            if(Cutpoint.size()>0 && Cur_range.first-Cutpoint.back().second<=MIN_GAP_SPAN)//Merge
                Cutpoint.back().second = Cur_range.second;
            else//NewRange
                Cutpoint.emplace_back(Cur_range);
        }
        return;
    }
    void 
    ReadGroup::GetCvgStats( std::map<uint32_t,int32_t>&coverage,
                            const std::map<uint32_t,uint32_t>& Boundary,
                            std::map<uint32_t,std::vector<double>>& CvgStats)
    {
        //Merging Boundary into coverage
        for(auto bound_iter=Boundary.begin();bound_iter!=Boundary.end();bound_iter++)
        {
            if(coverage.count(bound_iter->first)) // Exist
                continue;
            else // Add to coverage
            {
                auto insert {coverage.emplace(bound_iter->first,0)};
                if(insert.first->first!=(coverage.begin())->first) // Add to coverage
                    insert.first->second = std::prev(insert.first)->second; //inherit the coverage from prev
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
                MaxCvg    = (Inneriter->second>MaxCvg)?Inneriter->second:MaxCvg;
                Zerofrac += (Inneriter->second==0)?(Inneriter2->first-Inneriter->first):0;
                MeanCvg  += (Inneriter2->first-Inneriter->first)*Inneriter->second;
                Inneriter++;
                Inneriter2++;
                if(Inneriter2==coverage.end())
                    break;
            }
            CvgStats[Bound_iter->first]=std::move(std::vector<double>{MaxCvg,
                                                                      StartCvg,
                                                                      Zerofrac/(next_cvg_iter->first-Cvg_iter->first),
                                                                      MeanCvg/(next_cvg_iter->first-Cvg_iter->first)});
            Bound_iter++;
            next_Bound_iter++;
        }
        return;
    }
    /*
     * Give each read a Type and calculate Type statistics
     */
    void
    ReadGroup::CalculateType()
    {
        ExonCoverage.assign(ExonBoundary.size(),0);
        for(auto read_index=0;read_index<ReadStart.size();read_index++)
        {
            //Return type of current read
            auto ReadType {std::move(GetType(ReadStart[read_index],ReadEnd[read_index],ExonCoverage,MIN_OVERLAP))};
            auto TypeIter {std::find(SGTypes.begin(),SGTypes.end(),ReadType)};

            if(TypeIter==SGTypes.end())//New Type
            {
                Read2Type.emplace_back(SGTypes.size());
                SGTypes.emplace_back(ReadType);                
                TypeCount.emplace_back(1);
                TypeDirection.emplace_back(Direction[read_index]);
                ValidType.emplace_back(ReadType.size()>0);
            }
            else //Existing type
            {
                auto TypeIndex {std::distance(SGTypes.begin(),TypeIter)};
                TypeCount[TypeIndex]++;
                TypeDirection[TypeIndex] += Direction[read_index];
                Read2Type.emplace_back(TypeIndex);
            }
        }
        return;
    }

    /* 
     * Assign current read to exons by overlapping proportion.
     * By doing so, current read is now assigned to a type (a.k.a an exon combination)
     */  
    std::vector<uint32_t>
    ReadGroup::GetType(const std::vector<uint32_t>& SegStart,
                       const std::vector<uint32_t>& SegEnd,
                       std::vector<uint32_t>& ExonCoverage,
                       const uint32_t& MIN_OVERLAP)
    {
        std::vector<uint32_t> ReadType;
        auto leftmost_exon_index {0};
        for(auto seg_index=0;seg_index<SegStart.size();seg_index++)
        {
            //Get ovelapping exons for current read (segment)
            for(auto exon_index=leftmost_exon_index;exon_index<ExonBoundary.size();exon_index++)
            {
                uint32_t OverlapDist {utils::GetOverLapping(SegStart[seg_index],
                                                            SegEnd[seg_index],
                                                            ExonBoundary[exon_index].first,
                                                            ExonBoundary[exon_index].second)};                
                if(OverlapDist>0)//Is overlapping
                {
                    if(OverlapDist<MIN_OVERLAP)//Not Large enough
                    {
                        //First read
                        if( seg_index==0 &&
                            SegEnd[seg_index]>ExonBoundary[exon_index].second &&
                            ExonBoundary[exon_index].second>SegStart[seg_index] &&
                            SegStart[seg_index]>ExonBoundary[exon_index].first)
                            continue;
                        //Last read
                        if( seg_index==SegStart.size()-1 &&
                            SegEnd[seg_index]<ExonBoundary[exon_index].second &&
                            ExonBoundary[exon_index].first<SegEnd[seg_index] &&
                            SegStart[seg_index]<ExonBoundary[exon_index].first)
                            continue;
                    }
                    ExonCoverage[exon_index]++;
                    ReadType.emplace_back(exon_index);
                    leftmost_exon_index = exon_index; // Cannot be exon++, since two segments might be mapped to same exon. 
                }
            }
        }
        return ReadType;
    }

    void
    ReadGroup::RemoveWeakExons(const double MIN_CVG_FRAC)
    {
        std::vector<bool> IncludedInJunction(ExonBoundary.size(),false);
        
        for(auto TypeIndex=0;TypeIndex<SGTypes.size();TypeIndex++)
        {
            if(SGTypes[TypeIndex].size()>1)
            {
                if(SGTypes[TypeIndex][0]!=SGTypes[TypeIndex][1]-1)
                    IncludedInJunction[SGTypes[TypeIndex][0]] = true;
                if(SGTypes[TypeIndex].back()!=SGTypes[TypeIndex][SGTypes[TypeIndex].size()-2]+1)
                    IncludedInJunction[SGTypes[TypeIndex].back()] = true;
                if(SGTypes[TypeIndex].size()>2)
                {
                    for(auto exon_index=1;exon_index<SGTypes[TypeIndex].size()-1;exon_index++)
                        IncludedInJunction[SGTypes[TypeIndex][exon_index]] = true;
                }
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
            if(ExonBoundary.size()>1)
            {
                if(exon_index==0)
                    Neighbor_Cvg = CvgStats[ExonBoundary[1].first][3];
                else if(exon_index==ExonBoundary.size()-1)
                    Neighbor_Cvg = CvgStats[ExonBoundary[ExonBoundary.size()-2].first][3];
                else
                    Neighbor_Cvg = CvgStats[ExonBoundary[exon_index+1].first][3]>CvgStats[ExonBoundary[exon_index-1].first][3]?
                                CvgStats[ExonBoundary[exon_index+1].first][3]:
                                CvgStats[ExonBoundary[exon_index-1].first][3];
            }
            if(CvgStats[ExonBoundary[exon_index].first][3]<=Neighbor_Cvg*MIN_CVG_FRAC 
               && !IncludedInJunction[exon_index])
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

        for(const auto& bound:ExonBoundary)
        {
            Boundary[bound.first]++;
            Boundary[bound.second+1]++;
        }
        GetCvgStats(Coverage,Boundary,CvgStats);
        return;
    }

    void
    ReadGroup::WriteStatsToFile(std::ofstream& ofs,
                                const TwoDimVec<uint32_t>& CandidateIsf,
                                const std::vector<double>& ExpLv,
                                const std::vector<int64_t>& IsfDir)
    {
        ofs << "=================== BEGIN OF ReadGroup "<< RG_index+1<<"-"<<SubRG_index+1<<" ==================="<<std::endl;
        ofs << "[Boundary]\n" << ChrName<<":("<<CurrentRange.first<<"-"<<CurrentRange.second<<"), Orientation:"<<Orientation<<std::endl;       
        ofs << "[ReadLen]\n"  << ReadLen<<std::endl;
        ofs << "[Total valid reads]\n"<<validSize()<<std::endl;
        ofs << "[Exon Information]\n";
        ofs << std::left << std::setw(15) << "Exon_Index "
                         << std::setw(15) << "Exon_Start " 
                         << std::setw(15) << "Exon_End " 
                         << std::setw(15) << "Exon_Length " 
                         << std::setw(15) << "Exon_Coverage " 
                         << std::setw(15) << "Coverage Stats "<< std::endl;

        for(auto exon_index=0;exon_index<ExonBoundary.size();exon_index++)
        {
            if(ValidExons[exon_index])
            {
                ofs<<std::left << std::setw(15) << exon_index << " "
                               << std::setw(15) << ExonBoundary[exon_index].first <<" "
                               << std::setw(15) << ExonBoundary[exon_index].second<<" "
                               << std::setw(15) << (ExonBoundary[exon_index].second-ExonBoundary[exon_index].first+1)<<" "
                               << std::setw(15) << ExonCoverage[exon_index]<<" "<< std::setw(15);

                auto Cvg_iter {CvgStats.find(ExonBoundary[exon_index].first)};
                for(const auto& stats:Cvg_iter->second)
                    ofs<<stats<<",";
                ofs<<std::endl;
            
            }
        }

        ofs <<"[SGTypes]"<< std::endl;
        ofs << std::left << std::setw(15) << "Index "
                         << std::setw(15) << "Type " 
                         << std::setw(15) << "TypeCount"
                         << std::setw(15) << "Direction " 
                         << std::endl;

        for(auto TypeIndex=0;TypeIndex<SGTypes.size();TypeIndex++)
        {
            if(!ValidType[TypeIndex])
                continue;
            ofs<<std::left << std::setw(15)<<TypeIndex;
            auto Exon_iter {SGTypes[TypeIndex].begin()};

            for(auto exon_index=0;exon_index<ExonBoundary.size();exon_index++)
            {
                if(ValidExons[exon_index])
                {
                    if(Exon_iter==SGTypes[TypeIndex].end())
                        ofs<<"0";
                    else if(*Exon_iter==exon_index)
                    {
                        ofs<<"1";
                        Exon_iter++;
                    }
                    else
                        ofs<<"0";
                }
            }
            
            ofs<< std::setw(15)<<" "<<TypeCount[TypeIndex];
            if(TypeDirection[TypeIndex]>0)
                ofs<<" +"<<std::endl;
            else if(TypeDirection[TypeIndex]<0)
                ofs<<" -"<<std::endl;
            else
                ofs<<" ."<<std::endl;
        }

        uint32_t nPETypes {0};
        std::map<range_type,std::vector<uint32_t>> PETypes_map;

        for(auto read_index=0;read_index<ReadStart.size();read_index++)
        {
            if(PairendTable[read_index]!=-1 && read_index<PairendTable[read_index])
            {
                uint32_t LeftType  {Read2Type[read_index]};
                uint32_t RightType {Read2Type[PairendTable[read_index]]};

                if(!(ValidType[LeftType]&&ValidType[RightType]))
                    continue;

                nPETypes++;
                PETypes_map[range_type(LeftType,RightType)].emplace_back(ReadStart[PairendTable[read_index]].front()-ReadEnd[read_index].back());
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

        ofs<<"[PETypes]\n"<<nPETypes<<std::endl;
        for(auto PE_Iter=PETypes_map.begin();PE_Iter!=PETypes_map.end();PE_Iter++)
        {
            ofs<<ValidSGIndex[PE_Iter->first.first]<<" "<<
                 ValidSGIndex[PE_Iter->first.second]<<" "<<
                 PE_Iter->second.size()<<" ";
            ofs<<std::accumulate(PE_Iter->second.begin(),PE_Iter->second.end(),0.0)/PE_Iter->second.size()<<std::endl;
        }

        ofs<<"[Assembled Isoforms]"<<std::endl;
        ofs << std::left << std::setw(15) << "Index "
                         << std::setw(30) << "Isoform " 
                         << std::setw(15) << "Expression Level(RPKM) " 
                         << std::endl;
        for(auto Isf_index=0;Isf_index<CandidateIsf.size();++Isf_index)
        {
            ofs<<std::left << std::setw(15) << Isf_index <<"\t";

            auto iter {CandidateIsf[Isf_index].begin()};
            for(auto exon_index=0;exon_index<ExonBoundary.size();++exon_index)
            {
                if(ValidExons[exon_index])
                {
                    if(*iter==exon_index)
                    {
                        ofs<<"1";
                        ++iter;
                    }
                    else
                        ofs<<"0";
                }
            }
            ofs<<std::left << std::setw(15) << "\t" << ExpLv[Isf_index]
                           << std::endl;
        }
        ofs << "=================== END OF ReadGroup "<< RG_index+1<<"-"<<SubRG_index+1<<" ==================="<<std::endl;

        return;
    }   


    void
    ReadGroup::WritePredToGTF(std::ofstream& ofs,
                              const TwoDimVec<uint32_t>& CandidateIsf,
                              const std::vector<double>& ExpLv,
                              const std::vector<int64_t>& IsfDir,
                              const std::vector<double>& IsoformProb)
    {
        auto Isf_cnt {0};
        auto dir {'.'};
        for(const auto& Isf:CandidateIsf)
        {
            //Ignore low expressed isoform
            if(ExpLv[Isf_cnt]<ISF_EXPLV_THRESHOLD||IsoformProb[Isf_cnt]<IsoLasso::Algorithm::EMConfig::MinProbCut)
            {
                Isf_cnt++;
                continue;
            }
            if(IsfDir[Isf_cnt]>0)
                dir='+';
            else if(IsfDir[Isf_cnt]<0)
                dir='-';
            else
                dir='.';
            ofs<<ChrName<<"\t"
               <<"IsoLasso_"<<version<<"\t"
               <<"transcript"<<"\t"
               <<ExonBoundary[Isf.front()].first<<"\t"
               <<ExonBoundary[Isf.back()].second<<"\t"
               <<"."<<"\t"
               <<dir<<"\t"
               <<"."<<"\t"
               <<"gene_id \""<<"ReadGroup "<<RG_index+1<<"-"<<SubRG_index+1<<"\"; "
               <<"transcript_id \""<<"Isf"<<RG_index+1<<"_"<<SubRG_index+1<<"_"<<Isf_cnt+1<<"\"; "
               <<"RPKM \""<<ExpLv[Isf_cnt]<<"\"; "
               <<std::endl;
            
            //Merge contiguous exons
            auto exon_iter {Isf.begin()};
            auto exon_start {ExonBoundary[*exon_iter].first}, exon_end {ExonBoundary[*exon_iter].second}, exon_idx {*exon_iter};
            auto NewExonFlag {true};
            ++exon_iter;
            while(exon_iter!=Isf.end())
            {
                if(ExonBoundary[*exon_iter].first==(exon_end+1)) // Contiguous exon
                {
                    exon_end = ExonBoundary[*exon_iter].second;
                    ++exon_iter;
                    NewExonFlag = false;
                }
                else
                {
                    ofs <<ChrName<<"\t"
                        <<"IsoLasso_"<<version<<"\t"
                        <<"exon"<<"\t"
                        <<exon_start<<"\t"
                        <<exon_end<<"\t"
                        <<"."<<"\t"
                        <<dir<<"\t"
                        <<"."<<"\t"
                        <<"gene_id \""<<"ReadGroup "<<RG_index+1<<"-"<<SubRG_index+1<<"\"; "
                        <<"transcript_id \""<<"Isf" <<RG_index+1<<"_"<<SubRG_index+1<<"_"<<Isf_cnt+1<<"\"; "
                        <<"exon_number "<<exon_idx+1<<"; "
                        <<std::endl;   

                    exon_start  = ExonBoundary[*exon_iter].first;
                    exon_end    = ExonBoundary[*exon_iter].second;
                    exon_idx    = *exon_iter;
                    NewExonFlag = true;
                        
                    ++exon_iter;            
                }
            }

            ofs<<ChrName<<"\t"
               <<"IsoLasso_"<<version<<"\t"
               <<"exon"<<"\t"
               <<exon_start<<"\t"
               <<exon_end<<"\t"
               <<"."<<"\t"
               <<dir<<"\t"
               <<"."<<"\t"
               <<"gene_id \""<<"ReadGroup "<<RG_index+1<<"-"<<SubRG_index+1<<"\"; "
               <<"transcript_id \""<<"Isf"<<RG_index+1<<"_"<<SubRG_index+1<<"_"<<Isf_cnt+1<<"\"; "
               <<"exon_number "<<exon_idx+1<<"; "
               <<std::endl;
            

            ++Isf_cnt;
        }
        return;
    }

    void
    ReadGroup::PostProcess()
    {
        std::vector<uint32_t> Invalid_Cnt(ExonBoundary.size(),0);
        for(auto exon_index=1;exon_index<ExonBoundary.size();exon_index++)
            Invalid_Cnt[exon_index]=Invalid_Cnt[exon_index-1]+(ValidExons[exon_index-1]==false);
                
        uint32_t ValidTypeidx {0};
        std::vector<uint32_t> ValidSGIndex(ValidType.size(),0);
        for(auto Type_index=0;Type_index<ValidType.size();Type_index++)
        {
            if(ValidType[Type_index])
                ValidSGIndex[Type_index] = ValidTypeidx++;
            else
                ValidSGIndex[Type_index] = -1;
        }

        for(auto TypeIndex=0;TypeIndex<SGTypes.size();)
        {
            if(ValidType[TypeIndex])
            {
                for(auto exon_index=0;exon_index<SGTypes[TypeIndex].size();exon_index++)
                    SGTypes[TypeIndex][exon_index]-=Invalid_Cnt[SGTypes[TypeIndex][exon_index]];
                
                TypeIndex++;
            }
            else
            {
                ValidType.erase(ValidType.begin()+TypeIndex);
                SGTypes.erase(SGTypes.begin()+TypeIndex);
                TypeCount.erase(TypeCount.begin()+TypeIndex);
                TypeDirection.erase(TypeDirection.begin()+TypeIndex);
            }
        }

        //PETypes
        std::map<range_type,std::vector<uint32_t>> PETypes_map;

        for(auto read_index=0;read_index<ReadStart.size();read_index++)
        {
            if(PairendTable[read_index]!=-1 && read_index<PairendTable[read_index])
            {
                uint32_t LeftType  {ValidSGIndex[Read2Type[read_index]]};
                uint32_t RightType {ValidSGIndex[Read2Type[PairendTable[read_index]]]};

                if(LeftType==-1||RightType==-1)
                    continue;

                PEReadCount++;
                PETypes_map[range_type(LeftType,RightType)].emplace_back(ReadStart[PairendTable[read_index]].front()-ReadEnd[read_index].back());
            }
        }

        if(PETypes_map.size()>0)
        {
            PETypes.assign(PETypes_map.size(),std::vector<uint32_t>());
            uint32_t PE_Index {0};
            for(auto PE_Iter=PETypes_map.begin();PE_Iter!=PETypes_map.end();PE_Iter++)
            {
                PETypes[PE_Index].emplace_back(ValidSGIndex[PE_Iter->first.first]);
                PETypes[PE_Index].emplace_back(ValidSGIndex[PE_Iter->first.second]);
            }
        }

        for(auto exon_index=0;exon_index<ExonBoundary.size();)
        {
            if(ValidExons[exon_index])
            {
                ExonStats.emplace_back(CvgStats.find(ExonBoundary[exon_index].first)->second);
                exon_index++;
            }
            else
            {
                ExonBoundary.erase(ExonBoundary.begin()+exon_index);
                ValidExons.erase(ValidExons.begin()+exon_index);
            }
        }
        ReadCount = std::accumulate(ValidRead.begin(),ValidRead.end(),0,[](uint32_t a,bool b){return a+(b==true);}); 

        fillRSMatrix();

        return;
    }


}//end of IsoLasso::format
namespace IsoLasso::utils
{
    void
    ProcessReadGroup(IsoLasso::format::ReadGroup RG)
    {

#ifdef  DEBUG
        auto Start_time {std::chrono::steady_clock::now()};
#endif
        //If segments within same read has too high coverage, remove them.
        RG.RemoveLongSpanReads(MAX_EXON_SPAN,MAX_JUNCTION_COVERAGE);
        //Split RG into SubRGs
        std::vector<format::ReadGroup> SubRGs;
        RG.SplitbyRangeSet(SubRGs,MIN_GAP_SPAN);

#ifdef  DEBUG
        IsoLasso::utils::ShowRunningTime(Start_time,"I.","Split By RangeSet");
#endif

        //Process each sub-readgroup
        for(auto& SubRG:SubRGs)
        {
#ifdef  DEBUG
            Start_time = std::chrono::steady_clock::now();
#endif 
            //Enumerate all boundaries
            SubRG.CalculateBound(MIN_JUNC_COV,MIN_GAP_SPAN);
            //Calculate SGType (Combinations of high quality exon combinations)
            SubRG.CalculateType();
            //Remove Exons with too small coverage and not included in any junction type
            SubRG.RemoveWeakExons(MIN_CVG_FRAC);
            //Remove Types with Invalid exons and those reads with invalid types and recalculate coverage statistics
            SubRG.CalculateValidExons();

            if(SubRG.validSize() < MIN_RG_SIZE)
                continue;
            SubRG.PostProcess();
#ifdef  DEBUG
            IsoLasso::utils::ShowRunningTime(Start_time,"II.","Pre-Processing");
#endif
            /* E-M Module */
            Algorithm::PredExpLevel(SubRG);
        }
        return;
    }


    

}//end of IsoLasso::utils