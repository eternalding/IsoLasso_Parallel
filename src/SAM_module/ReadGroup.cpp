#include <SAM_module/ReadGroup.hpp>
#include <utils/Commontype.hpp>

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
    SplitbyRangeSet(std::vector<IsoLasso::format::ReadGroup>& SubRGs,const uint32_t MIN_GAP_SPAN)
    {
        


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
    }

}//end of IsoLasso::utils
