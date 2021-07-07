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




        return;
    }


}// end of namespace IsoLasso::algorithm