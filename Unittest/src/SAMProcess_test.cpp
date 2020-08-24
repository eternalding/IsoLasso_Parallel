#include <gtest/gtest.h>
#include "SAM_module/SAMprocess.hpp"
#include "SAM_module/SAMrecord.hpp"
#include "SAM_module/ReadGroup.hpp"
#include <string>
#include <sstream>
#include <fstream>

using namespace IsoLasso::format;
using namespace IsoLasso::utils;

/*
TEST(SamProcess, SplitReadGroup)
{
    IsoLasso::format::Header_record HRecord;
    IsoLasso::format::Sam_record Record;
    IsoLasso::format::ReadGroup RG;
    range_type CurrentRange(0,0);
    std::ostringstream os;

    std::ifstream fin("../Unittest/TestData/example2_top1000000.sam",std::ios::in);

    EXPECT_EQ(fin.is_open(),true);

    while(fin>>HRecord);
    while(fin>>Record)
    {
        if(Record.Pos<CurrentRange.first)
            std::__throw_invalid_argument("Using unsorted SAM file! Aborting...");

        if(Record.ValidBit)
        {
            //Add current record to current ReadGroup
            RG.AddRecord(Record,CurrentRange);
        }
    }

    EXPECT_EQ(RG.PairendTable.size(),27);

    std::vector<long> expected_result(27,-1);
    expected_result[0]  = 1;
    expected_result[1]  = 0;

    expected_result[6]  = 8;
    expected_result[8]  = 6;

    expected_result[7]  = 13;
    expected_result[13] = 7;

    expected_result[16] = 17;
    expected_result[17] = 16;

    expected_result[20] = 26;
    expected_result[26] = 20;

    for(auto i=0;i<expected_result.size();i++)
        EXPECT_EQ(RG.PairendTable[i],expected_result[i]);
};
*/