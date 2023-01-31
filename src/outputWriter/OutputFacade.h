/*
 * OutputFacade.h
 *
 * Created: 01.11.2022
 * Author:  borisov
 */

#pragma once

#include "./XYZWriter.h"
#include "./VTKWriter.h"
#include "./CheckpointWriter.h"
#include "./CSVWriter.h"
#include "../model/ProgramParameters.h"
#include "spdlog/spdlog.h"

/**
 * @brief: Facade which handles all the output functionalities
 */
class OutputFacade
{
private:
    ProgramParameters* programParameters;
    std::string prefix; 
    outputWriter::XYZWriter xyzWriter;
    outputWriter::VTKWriter vtkWriter;
    outputWriter::CheckpointWriter checkpointWriter;
    outputWriter::CSVWriter csvWriter = outputWriter::CSVWriter("../csv");

    /**
     * A spdlog logger, which logs the logic of the program flow
     */
    std::shared_ptr<spdlog::logger> _logicLogger;

    /**
     * a speedlog logger which logs construction and destruction of particles
     */
    std::shared_ptr<spdlog::logger> _memoryLogger;

    /**
     * @brief creates a directory
     * @param path of the directory to be created
     */
    void createDirectory(std::string path);

    /**
     * @brief deletes a directory
     * @param path of the directory to be deleted
     */
    void removeDirectory(std::string path);

public:
    /**
     * @brief Constructs a new Output Facade object and creates folders which are needed for output
     * @param programParameters program parameters
     */
    OutputFacade(ProgramParameters* programParameters);

    ~OutputFacade();

    /**
     * @brief writes .xyz files, containing the simulation's data, to the build/outputXYZ/ folder
     * @param iteration number of the iteration, which is used to generate a unique filename
     */
    void outputXYZ(int iteration);

    /**
     * @brief writes .vtu files, containing the simulation's data, to the build/outputVTK/ folder
     * @param iteration number of the iteration, which is used to generate a unique filename
     */
    void outputVTK(int iteration);

    /**
     * @brief creates a checkpoint with the state of all particles
     */
    void createCheckpoint();

    /**
     * @brief writes the data to csv file (append)
     * 
     * @param data the data to write
     * @param avg the total sum of all particles ever measured in a bin
     */
    void writeCSV(std::vector<int> data, std::vector<int> avg);
};