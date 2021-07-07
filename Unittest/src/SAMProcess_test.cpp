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

inline IsoLasso::format::Args
RunBinary(char** argv,const int argc)
{
    IsoLasso::format::Args arguments;
    IsoLasso::utils::ParseArg(argc,argv,arguments);

    return arguments;
}


TEST(SamProcess, ReadGroupCount)
{

    char* argv[] = {(char*)"./IsoLasso_Parallel",
                    (char*)"-i",
                    (char*)"../Unittest/TestData/example2_top10000.sam",
                    (char*)"-o",
                    (char*)"output"};
    int argc = 5;

    IsoLasso::format::Args arguments = RunBinary(argv,argc);
    auto [RGcount,readcount]         = IsoLasso::utils::ReadSamFile(arguments);

    EXPECT_EQ(RGcount,2);
    EXPECT_EQ(readcount,9901);

};

TEST(SamProcess, HeaderAtTheEnd)
{
    char* argv[] = {(char*)"./IsoLasso_Parallel",
                    (char*)"-i",
                    (char*)"../Unittest/TestData/3exon_small_40_60_depth20.sam",
                    (char*)"-o",
                    (char*)"output"};
    int argc = 5;
    IsoLasso::format::Args arguments=RunBinary(argv,argc);

    auto [RGcount,readcount] = IsoLasso::utils::ReadSamFile(arguments);

    EXPECT_EQ(RGcount,2);
    EXPECT_EQ(readcount,2013);

};

TEST(SamProcess, SplitBySubGroup)
{
    char* argv[] = {(char*)"./IsoLasso_Parallel",
                    (char*)"-i",
                    (char*)"../Unittest/TestData/example2_top10000.sam",
                    (char*)"-o",
                    (char*)"output"};
    int argc = 5;

    IsoLasso::format::Args arguments=RunBinary(argv,argc);
    
    auto [RGcount,readcount] = IsoLasso::utils::ReadSamFile(arguments);

    EXPECT_EQ(RGcount,2);
    EXPECT_EQ(readcount,9901);
    
};

TEST(SamProcess, ExonBoundary)
{
    IsoLasso::format::ReadGroup RG;
    RG.ReadStart.emplace_back(std::initializer_list<uint32_t>{10,60,140,270,580,651,660,670,680}); 
    RG.ReadEnd.emplace_back(std::initializer_list<uint32_t>{100,120,180,300,650,700,700,700,700}); 
    RG.CurrentRange=RG.getRange();
    RG.ValidRead.assign(RG.ReadStart.size(),true);
    RG.ReadLen = 75;

    RG.CalculateBound(MIN_JUNC_COV,MIN_GAP_SPAN);

    EXPECT_EQ(RG.CurrentRange.first,10);
    EXPECT_EQ(RG.CurrentRange.second,700);

    std::vector<range_type> GT_range {std::pair<uint32_t,uint32_t>(10,59),
                                      std::pair<uint32_t,uint32_t>(60,100),
                                      std::pair<uint32_t,uint32_t>(101,120),
                                      std::pair<uint32_t,uint32_t>(121,139),
                                      std::pair<uint32_t,uint32_t>(140,180),
                                      std::pair<uint32_t,uint32_t>(181,269),
                                      std::pair<uint32_t,uint32_t>(270,300),
                                      std::pair<uint32_t,uint32_t>(301,579),
                                      std::pair<uint32_t,uint32_t>(580,650),
                                      std::pair<uint32_t,uint32_t>(651,659),
                                      std::pair<uint32_t,uint32_t>(660,669),
                                      std::pair<uint32_t,uint32_t>(670,679),
                                      std::pair<uint32_t,uint32_t>(680,700)};

    EXPECT_EQ(RG.ExonBoundary,GT_range);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[0].first][0],1);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[1].first][0],2);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[2].first][0],1);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[3].first][0],0);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[4].first][0],1);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[5].first][0],0);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[6].first][0],1);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[7].first][0],0);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[8].first][0],1);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[9].first][0],1);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[10].first][0],2);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[11].first][0],3);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[12].first][0],4);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[0].first][3],1);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[1].first][3],2);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[2].first][3],1);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[3].first][3],0);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[4].first][3],1);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[5].first][3],0);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[6].first][3],1);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[7].first][3],0);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[8].first][3],1);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[9].first][3],1);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[10].first][3],2);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[11].first][3],3);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[12].first][3],4);
}

TEST(SamProcess, SGTypes)
{
    IsoLasso::format::ReadGroup RG;
    RG.ReadStart.emplace_back(std::initializer_list<uint32_t>{10,60,140,270,580,651,660,670,680}); 
    RG.ReadEnd.emplace_back(std::initializer_list<uint32_t>{100,120,180,300,650,700,700,700,700}); 
    RG.Direction.assign(RG.ReadStart.size(),0);

    RG.CurrentRange=RG.getRange();
    RG.ValidRead.assign(RG.ReadStart.size(),true);

    RG.CalculateBound(MIN_JUNC_COV,MIN_GAP_SPAN);
    RG.CalculateType();
    
    /*
    EXPECT_EQ(RG.CurrentRange.first,10);
    EXPECT_EQ(RG.CurrentRange.second,700);

    std::vector<range_type> GT_range {std::pair<uint32_t,uint32_t>(10,59),
                                      std::pair<uint32_t,uint32_t>(60,100),
                                      std::pair<uint32_t,uint32_t>(101,120),
                                      std::pair<uint32_t,uint32_t>(121,139),
                                      std::pair<uint32_t,uint32_t>(140,141),
                                      std::pair<uint32_t,uint32_t>(142,180),
                                      std::pair<uint32_t,uint32_t>(181,269),
                                      std::pair<uint32_t,uint32_t>(270,271),
                                      std::pair<uint32_t,uint32_t>(272,300),
                                      std::pair<uint32_t,uint32_t>(301,579),
                                      std::pair<uint32_t,uint32_t>(580,650),
                                      std::pair<uint32_t,uint32_t>(651,659),
                                      std::pair<uint32_t,uint32_t>(660,669),
                                      std::pair<uint32_t,uint32_t>(670,679),
                                      std::pair<uint32_t,uint32_t>(680,700)};

    EXPECT_EQ(RG.ExonBoundary,GT_range);

    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[0].first][0],1);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[1].first][0],2);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[2].first][0],1);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[3].first][0],0);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[4].first][0],1);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[5].first][0],1);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[6].first][0],0);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[7].first][0],1);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[8].first][0],1);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[9].first][0],0);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[10].first][0],1);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[11].first][0],1);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[12].first][0],2);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[13].first][0],3);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[14].first][0],4);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[0].first][3],1);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[1].first][3],2);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[2].first][3],1);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[3].first][3],0);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[4].first][3],1);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[5].first][3],1);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[6].first][3],0);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[7].first][3],1);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[8].first][3],1);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[9].first][3],0);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[10].first][3],1);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[11].first][3],1);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[12].first][3],2);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[13].first][3],3);
    EXPECT_EQ(RG.CvgStats[RG.ExonBoundary[14].first][3],4);
    */
}