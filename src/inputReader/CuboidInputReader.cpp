/*
 * CuboidInputReader.cpp
 *
 *  Created on: 16.11.2022
 *      Author: borisov
 */

#include "./CuboidInputReader.h"
#include "../utils/ArrayUtils.h"
#include "../utils/ParticleGenerator.h"

#include <string>
#include <fstream>
#include <sstream>

void CuboidInputReader::readInput(ProgramParameters &programParameters, const char *filename)
{
    // Variables to read in
    std::array<double, 3> x;
    std::array<int, 3> n;
    double h;
    double m;
    std::array<double, 3> v;
    double epsilon;
    double sigma;
    int type;
    int fixed_int;
    bool fixed_bool;

    std::ifstream input_file(filename);
    std::string tmp_string;
    auto particleContainer = programParameters.getParticleContainer();

    if (input_file.is_open())
    {
        // skip comments plus cuboid file indicator
        while (tmp_string.empty() or tmp_string[0] == '#' or tmp_string[0] == '$')
        {
            getline(input_file, tmp_string);
            getLogicLogger()->info("Read line: {}", tmp_string);
        }
        // tmp_string now contains the x
        std::istringstream datastream(tmp_string);

        for (auto &xi : x)
        {
            datastream >> xi;
        }
        datastream.clear();

        // get next line which contains the n
        getline(input_file, tmp_string);
        getLogicLogger()->info("Read line: {}", tmp_string);
        datastream.str(tmp_string);

        for (auto &ni : n)
        {
            datastream >> ni;
        }
        datastream.clear();

        // get next line which contains the h
        getline(input_file, tmp_string);
        getLogicLogger()->info("Read line: {}", tmp_string);
        datastream.str(tmp_string);

        datastream >> h;
        datastream.clear();

        // get next line which contains the m
        getline(input_file, tmp_string);
        getLogicLogger()->info("Read line: {}", tmp_string);
        datastream.str(tmp_string);

        datastream >> m;
        datastream.clear();

        // get next line which contains the v
        getline(input_file, tmp_string);
        getLogicLogger()->info("Read line: {}", tmp_string);
        datastream.str(tmp_string);

        for (auto &vi : v)
        {
            datastream >> vi;
        }
        datastream.clear();

        // get next line which contains the epsilon
        getline(input_file, tmp_string);
        getLogicLogger()->info("Read line: {}", tmp_string);
        datastream.str(tmp_string);

        datastream >> epsilon;
        datastream.clear();

        // get next line which contains the sigma
        getline(input_file, tmp_string);
        getLogicLogger()->info("Read line: {}", tmp_string);
        datastream.str(tmp_string);

        datastream >> sigma;
        datastream.clear();

        // get next line wich contains the type
        getline(input_file, tmp_string);
        getLogicLogger()->info("Read line: {}", tmp_string);
        datastream.str(tmp_string);

        datastream >> type;
        datastream.clear();

        // get next line wich contains the fixed bool
        getline(input_file, tmp_string);
        getLogicLogger()->info("Read line: {}", tmp_string);
        datastream.str(tmp_string);

        datastream >> fixed_int;
        datastream.clear();

        if (fixed_int == 1) {
            fixed_bool = true;
        } else {
            fixed_bool = false;
        }
    }
    else
    {
        getLogicLogger()->error("Error: could not open file {}", filename);
        exit(-1);
    }

    // the membrane simulation only works through xml
    // so hardcoded stiffness and averageBondlength
    std::unique_ptr<Cuboid> cuboid = std::make_unique<Cuboid>(Cuboid(x, n, h, m, v, epsilon, sigma, type, 1, 1, fixed_bool));
    ParticleGenerator::generateCuboid(*particleContainer.get(), *cuboid, programParameters.getMembrane());

}