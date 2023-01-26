/*
 *  CSVWriter.cpp
 *
 *  Created on: 26.01.2023
 *      Author: borisov
 */

#include "./CSVWriter.h"
#include "../utils/File.h"

#include <iostream>
#include <fstream>

namespace outputWriter {
    // Constructor
    CSVWriter::CSVWriter(std::string output_path) {
        _logicLogger = spdlog::get("output_logger");
        _memoryLogger = spdlog::get("memory_logger");

        this->output_path = output_path + "/MD_DensityVelocity__" + File::getDateTime() + ".csv";
        this->count = 1;

        _memoryLogger->info("CSVWriter generated!");
    }

    CSVWriter::~CSVWriter() {
        _memoryLogger->info("CSVWriter destructed!");
    }

    void CSVWriter::writeData(std::vector<int> data, std::vector<int> avg) {
        std::ofstream file;
        file.open(output_path, std::ios::app);

        if (count == 1) {
            file << "Bins,";
            for (long unsigned int i = 0; i < data.size(); i++) {
                file << "Bin " << i << ",";
            }
            file << std::endl;
        }

        file << "Particles per bin:,";

        for (int i: data) {
            file << i << ",";
        }
        file << std::endl;

        file << "Average particles per bin:,";

        for (int i: avg) {
            file << (float)i / count << ",";
        }
        file << std::endl;

        count++;
        file.close();
    }
}