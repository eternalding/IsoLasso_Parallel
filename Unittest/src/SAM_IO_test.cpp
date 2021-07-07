#include <gtest/gtest.h>
#include "SAM_module/SAMprocess.hpp"
#include "SAM_module/SAMrecord.hpp"
#include "SAM_module/ReadGroup.hpp"
#include "EM_module/PredExpLevel.hpp"
#include "EM_module/EMalgorithm.hpp"
#include <string>
#include <sstream>
#include <fstream>

using namespace IsoLasso::format;
using namespace IsoLasso::utils;
using namespace IsoLasso::Algorithm;

TEST(Gtest, Gtest_initialization) 
{
  EXPECT_TRUE(true) << "This should not fail";
};


TEST(Arugment_test, ArgParsing)
{
    char* argv[] = {(char*)"./IsoLasso_Parallel",(char*)"-i",(char*)"test.sam",(char*)"-o",(char*)"output"};
    int argc = 5;

    IsoLasso::format::Args arguments;

    IsoLasso::utils::ParseArg(argc,argv,arguments);

    EXPECT_EQ(arguments.SAMFILE,"test.sam");
    EXPECT_EQ(arguments.OUTPUT_PREFIX,"output");

};

TEST(SAM_IO_test, Check_all_fields){

    IsoLasso::format::Sam_record Record;

    std::istringstream iss{R"(K00208:8901006:YAP006:8:2212:14397:16489        163     chr1    3044937    60       101M    =       3044987 151     GTGCCTAGAGGCTGCCTGGGGCTGAGAAAAGAGAAAAACAAACCTGGGTATGCCTCGTAGTTAAAACATTCCTGGGAACATCTTGACCATAAGATAAAGGG   AAFFFJJJFJJJJJJJ<JJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJFJJJJAJJJJJFFJJJJJJJJJJJJJJ<FFJFJJJJJJJJAJJFJJJJAJFJ   NH:i:6  HI:i:1  nM:i:1  AS:i:198        RG:Z:YAP006L8_TGACCAA XS:A:+)"};

    iss>>Record;

    EXPECT_EQ(Record.QName,"K00208:8901006:YAP006:8:2212:14397:16489");
    EXPECT_EQ(Record.Flag,163);
    EXPECT_EQ(Record.RName,"chr1");
    EXPECT_EQ(Record.Pos,3044937);
    EXPECT_EQ(Record.MapQ,60);
    EXPECT_EQ(Record.CIGAR,"101M");
    EXPECT_EQ(Record.RNext,"=");
    EXPECT_EQ(Record.PNext,3044987);
    EXPECT_EQ(Record.TLen,151);
    EXPECT_EQ(Record.Seq,"GTGCCTAGAGGCTGCCTGGGGCTGAGAAAAGAGAAAAACAAACCTGGGTATGCCTCGTAGTTAAAACATTCCTGGGAACATCTTGACCATAAGATAAAGGG");
    EXPECT_EQ(Record.Qual,"AAFFFJJJFJJJJJJJ<JJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJFJJJJAJJJJJFFJJJJJJJJJJJJJJ<FFJFJJJJJJJJAJJFJJJJAJFJ");
    EXPECT_EQ(Record.SpliceDir,1);

};

TEST(SAM_IO_test, Print_one_record){

    IsoLasso::format::Sam_record Record;

    std::istringstream iss{R"(K00208:8901006:YAP006:8:2212:14397:16489        163     chr1    3044937    60       101M    =       3044987 151     GTGCCTAGAGGCTGCCTGGGGCTGAGAAAAGAGAAAAACAAACCTGGGTATGCCTCGTAGTTAAAACATTCCTGGGAACATCTTGACCATAAGATAAAGGG   AAFFFJJJFJJJJJJJ<JJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJFJJJJAJJJJJFFJJJJJJJJJJJJJJ<FFJFJJJJJJJJAJJFJJJJAJFJ   NH:i:6  HI:i:1  nM:i:1  AS:i:198        RG:Z:YAP006L8_TGACCAA XS:A:+)"};

    iss>>Record;

    std::ostringstream os;

    os<<Record;

    EXPECT_EQ(os.str(),"K00208:8901006:YAP006:8:2212:14397:16489\t163\tchr1\t3044937\t60\t101M\t=\t3044987\t151\tGTGCCTAGAGGCTGCCTGGGGCTGAGAAAAGAGAAAAACAAACCTGGGTATGCCTCGTAGTTAAAACATTCCTGGGAACATCTTGACCATAAGATAAAGGG\tAAFFFJJJFJJJJJJJ<JJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJFJJJJAJJJJJFFJJJJJJJJJJJJJJ<FFJFJJJJJJJJAJJFJJJJAJFJ\t1");

};

TEST(SAM_IO_test, Print_Multiple_record){

    IsoLasso::format::Sam_record Record;

    std::istringstream iss{R"(r001    163 ref 7   30  8M4I4M1D3M  =   37  39  TTAGATAAAGAGGATACTG *   XX:B:S,12561,2,20,112
r002    0   ref 9   30  1S2I6M1P1I1P1I4M2I  *   0   0   AAAAGATAAGGGATAAA   *
r003    0   ref 9   30  5H6M    *   0   0   AGCTAA  *
r004    0   ref 16  30  6M14N1I5M   *   0   0   ATAGCTCTCAGC    *
r003    16  ref 29  30  6H5M    *   0   0   TAGGC   *
r001    83  ref 37  30  9M  =   7   -39 CAGCGCCAT   *)"};

    std::ostringstream os;

    while(iss>>Record)
      os<<Record;

    std::string Expected_result="r001\t163\tref\t7\t30\t8M4I4M1D3M\t=\t37\t39\tTTAGATAAAGAGGATACTG\t*\t0\
r002\t0\tref\t9\t30\t1S2I6M1P1I1P1I4M2I\t*\t0\t0\tAAAAGATAAGGGATAAA\t*\t0\
r003\t0\tref\t9\t30\t5H6M\t*\t0\t0\tAGCTAA\t*\t0\
r004\t0\tref\t16\t30\t6M14N1I5M\t*\t0\t0\tATAGCTCTCAGC\t*\t0\
r003\t16\tref\t29\t30\t6H5M\t*\t0\t0\tTAGGC\t*\t0\
r001\t83\tref\t37\t30\t9M\t=\t7\t-39\tCAGCGCCAT\t*\t0";

    EXPECT_EQ(os.str(),Expected_result);

};

TEST(SAM_IO_test, Skip_Header){

    IsoLasso::format::Header_record HRecord;
    IsoLasso::format::Sam_record Record;

    std::istringstream iss{R"(@HD VN:1.3  SO:coordinate
@SQ SN:ref  LN:45
@SQ SN:ref2 LN:40
r001    163 ref 7   30  8M4I4M1D3M  =   37  39  TTAGATAAAGAGGATACTG *   XX:B:S,12561,2,20,112
r002    0   ref 9   30  1S2I6M1P1I1P1I4M2I  *   0   0   AAAAGATAAGGGATAAA   *
r003    0   ref 9   30  5H6M    *   0   0   AGCTAA  *
r004    0   ref 16  30  6M14N1I5M   *   0   0   ATAGCTCTCAGC    *
r003    16  ref 29  30  6H5M    *   0   0   TAGGC   *
r001    83  ref 37  30  9M  =   7   -39 CAGCGCCAT   *)"};

    std::ostringstream os;

    while(iss>>HRecord);

    while(iss>>Record)
      os<<Record;

    std::string Expected_result="r001\t163\tref\t7\t30\t8M4I4M1D3M\t=\t37\t39\tTTAGATAAAGAGGATACTG\t*\t0\
r002\t0\tref\t9\t30\t1S2I6M1P1I1P1I4M2I\t*\t0\t0\tAAAAGATAAGGGATAAA\t*\t0\
r003\t0\tref\t9\t30\t5H6M\t*\t0\t0\tAGCTAA\t*\t0\
r004\t0\tref\t16\t30\t6M14N1I5M\t*\t0\t0\tATAGCTCTCAGC\t*\t0\
r003\t16\tref\t29\t30\t6H5M\t*\t0\t0\tTAGGC\t*\t0\
r001\t83\tref\t37\t30\t9M\t=\t7\t-39\tCAGCGCCAT\t*\t0";

    EXPECT_EQ(os.str(),Expected_result);

};

TEST(SAM_IO_test, Unmapped_record){

    IsoLasso::format::Header_record HRecord;
    IsoLasso::format::Sam_record Record;

    std::istringstream iss{R"(@HD VN:1.3  SO:coordinate
@SQ SN:ref  LN:45
@SQ SN:ref2 LN:40
r001    163 ref 7   30  8M4I4M1D3M  =   37  39  TTAGATAAAGAGGATACTG *   XX:B:S,12561,2,20,112
r002    0   ref 9   30  1S2I6M1P1I1P1I4M2I  *   0   0   AAAAGATAAGGGATAAA   *
r003    0   ref 9   30  5H6M    *   0   0   AGCTAA  *
r004    4   ref 16  30  6M14N1I5M   *   0   0   ATAGCTCTCAGC    *
r003    16  ref 29  30  6H5M    *   0   0   TAGGC   *
r001    83  ref 37  30  9M  =   7   -39 CAGCGCCAT   *)"};

    std::ostringstream os;

    while(iss>>HRecord);

    while(iss>>Record)
    {
      if(Record.ValidBit) 
        os<<Record;
    }

    std::string Expected_result="r001\t163\tref\t7\t30\t8M4I4M1D3M\t=\t37\t39\tTTAGATAAAGAGGATACTG\t*\t0\
r002\t0\tref\t9\t30\t1S2I6M1P1I1P1I4M2I\t*\t0\t0\tAAAAGATAAGGGATAAA\t*\t0\
r003\t0\tref\t9\t30\t5H6M\t*\t0\t0\tAGCTAA\t*\t0\
r003\t16\tref\t29\t30\t6H5M\t*\t0\t0\tTAGGC\t*\t0\
r001\t83\tref\t37\t30\t9M\t=\t7\t-39\tCAGCGCCAT\t*\t0";

    EXPECT_EQ(os.str(),Expected_result);

};

TEST(SAM_utils_test, CIGAR_parsing){

    IsoLasso::format::Sam_record Record;
std::istringstream iss{R"(r001    163 ref 7   30  8M4I4M1D3M  =   37  39  TTAGATAAAGAGGATACTG *   XX:B:S,12561,2,20,112
r002    0   ref 9   30  1S2I6M1P1I1P1I4M2I  *   0   0   AAAAGATAAGGGATAAA   *
r003    0   ref 9   30  5H6M    *   0   0   AGCTAA  *
r004    0   ref 16  30  6M14N1I5M   *   0   0   ATAGCTCTCAGC    *
r003    16  ref 29  30  6H5M    *   0   0   TAGGC   *
SRR959756.117    16    RNU6-1295P    110    9    19S25M    *    0    0    TAGAGACGGGGTCTCGCTATGTTGCTCAGGCTGGAGTGCAGTGG    VW[XY__c_f]^a^^a[dad[YYZ\V]]b]YZa\Y^[^^afaaa    AS:i:50    XS:i:42    XN:i:0    XM:i:0    XO:i:0    XG:i:0    NM:i:0    MD:Z:25    YT:Z:UU)"};

    std::ostringstream os;

    while(iss>>Record)
    {
      if(Record.ValidBit) 
      {
          for(auto start:Record.SegmentStart) 
            os<<start<<" ";
          os<<",";
          for(auto end:Record.SegmentEnd) 
            os<<end<<" ";
          os<<"\n";
      }
    }

    std::string expected_result{R"(7 ,22 
9 ,18 
9 ,14 
16 36 ,21 40 
29 ,33 
110 ,134 
)"};

    EXPECT_EQ(os.str(),expected_result);
}

TEST(SAM_utils_test, GetEfficientLength){

    IsoLasso::format::Sam_record Record;
std::istringstream iss{R"(r001    163 ref 7   30  8M4I4M1D3M  =   37  39  TTAGATAAAGAGGATACTG *   XX:B:S,12561,2,20,112
r002    0   ref 9   30  1S2I6M1P1I1P1I4M2I  *   0   0   AAAAGATAAGGGATAAA   *
r003    0   ref 9   30  5H6M    *   0   0   AGCTAA  *
r004    0   ref 16  30  6M14N1I5M   *   0   0   ATAGCTCTCAGC    *
r003    16  ref 29  30  6H5M    *   0   0   TAGGC   *
SRR959756.117    16    RNU6-1295P    110    9    19S25M    *    0    0    TAGAGACGGGGTCTCGCTATGTTGCTCAGGCTGGAGTGCAGTGG    VW[XY__c_f]^a^^a[dad[YYZ\V]]b]YZa\Y^[^^afaaa    AS:i:50    XS:i:42    XN:i:0    XM:i:0    XO:i:0    XG:i:0    NM:i:0    MD:Z:25    YT:Z:UU)"};

    std::ostringstream os;

    while(iss>>Record)
    {
      if(Record.ValidBit) 
      {
            os<<IsoLasso::utils::GetEfficientLen(Record)<<"\n";
      }
    }

    std::string expected_result{R"(16
10
6
11
5
25
)"};

    EXPECT_EQ(os.str(),expected_result);
}


TEST(SAM_utils_test, Paired_end_test){

    IsoLasso::format::Header_record HRecord;
    IsoLasso::format::Sam_record Record;
    IsoLasso::format::ReadGroup RG;
    range_type CurrentRange(0,0);
    std::ostringstream os;

    std::ifstream fin("../Unittest/TestData/example_3_paired.sam",std::ios::in);

    EXPECT_EQ(fin.is_open(),true);

    while(fin>>HRecord);
    while(fin>>Record)
    {
        if(Record.Pos<CurrentRange.first)
            std::__throw_invalid_argument("Using unsorted SAM file! Aborting...");

        if(Record.ValidBit)
        {
            RG.AddRecord(Record);
        }
    }

    EXPECT_EQ(RG.PairendTable.size(),27);

    std::vector<long> expected_result(27,-1);
    expected_result[0]  = 1;
    expected_result[1]  = 0;

    expected_result[2]  = 3;
    expected_result[3]  = 2;

    expected_result[4]  = 5;
    expected_result[5]  = 4;

    expected_result[6]  = 8;
    expected_result[8]  = 6;

    expected_result[7]  = 13;
    expected_result[13] = 7;

    expected_result[10] = 14;
    expected_result[14] = 10;

    expected_result[11] = 12;
    expected_result[12] = 11;

    expected_result[16] = 17;
    expected_result[17] = 16;

    expected_result[20] = 26;
    expected_result[26] = 20;

    for(auto i=0;i<expected_result.size();i++)
        EXPECT_EQ(RG.PairendTable[i],expected_result[i]);        
}
