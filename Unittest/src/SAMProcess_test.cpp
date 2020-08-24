#include <gtest/gtest.h>
#include "SAM_module/SAMprocess.hpp"
#include "SAM_module/SAMrecord.hpp"
#include "SAM_module/ReadGroup.hpp"
#include <string>
#include <sstream>
#include <fstream>

using namespace IsoLasso::format;
using namespace IsoLasso::utils;


TEST(SamProcess, SplitReadGroup)
{

    char* argv[] = {(char*)"./IsoLasso_Parallel",(char*)"-i",(char*)"testdata/example2_top10000.sam",(char*)"-o",(char*)"output"};
    int argc = 5;

    IsoLasso::format::Args arguments;

    IsoLasso::utils::ParseArg(argc,argv,arguments);

    auto [RGcount,readcount] = IsoLasso::utils::ReadSamFile(arguments);


    EXPECT_EQ(RGcount,7);
    EXPECT_EQ(readcount,9908);

};
