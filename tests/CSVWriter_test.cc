#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <fstream>
#include <cstdio>

#include "../src/outputWriter/CSVWriter.h"

// basic test for writing csv output
TEST(CSVWriter, WriteCsv)
{
    auto csvWriter = outputWriter::CSVWriter("../tests/csv");
    std::vector<int> data = {10, 10, 10, 10, 10};
    std::vector<int> avg = {20, 20, 20, 20, 20};

    csvWriter.writeData(data, avg);

    std::ifstream file(csvWriter.getOutputPath());
    std::string tmp_string;

    if (file.is_open()) {
        getline(file, tmp_string);
        EXPECT_EQ(tmp_string, "Bins,Bin 0,Bin 1,Bin 2,Bin 3,Bin 4,");
        getline(file, tmp_string);
        EXPECT_EQ(tmp_string, "Particles per bin:,10,10,10,10,10,");
        getline(file, tmp_string);
        EXPECT_EQ(tmp_string, "Average particles per bin:,20,20,20,20,20,");
    } else {
        FAIL();
    }
    std::remove(csvWriter.getOutputPath().c_str());
}