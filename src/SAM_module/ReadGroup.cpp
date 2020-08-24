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


}//end of IsoLasso::format

namespace IsoLasso::utils
{

    void
    ProcessReadGroup(IsoLasso::format::ReadGroup& RG,const range_type& current_range)
    {

    }

}//end of IsoLasso::utils
