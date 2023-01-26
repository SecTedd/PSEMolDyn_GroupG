/*
 *  CSVWriter.h
 *
 *  Created on: 26.01.2023
 *      Author: borisov
 */

#pragma once

#include "spdlog/spdlog.h"

namespace outputWriter {
    /**
     * @brief class that writes data to a csv file
     */
    class CSVWriter {
    private:
        // the output path
        std::string output_path; 

        // counter for how many times output was written (used in averages per bin calculation)
        int count; 

        // a speedlog logger which logs the logic flow of the simulation
        std::shared_ptr<spdlog::logger> _logicLogger;
        
        //a speedlog logger which logs construction and destruction of particles
        std::shared_ptr<spdlog::logger> _memoryLogger;

    public:
        /**
         * @brief Construct a new CSVWriter object
         * 
         * @param output_path the path for the output file
         */
        CSVWriter(std::string output_path);

        // destructor
        ~CSVWriter();

        /**
         * @brief writes the csv file to output_path
         * 
         * @param data the integer data to write in the file
         * @param avg the total sum of all particles ever measured in a bin
         */
        void writeData(std::vector<int> data, std::vector<int> avg);
    };
}