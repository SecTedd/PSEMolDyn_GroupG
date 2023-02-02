/*
 * SphereInputReader.cpp
 *
 *  Created on: 04.12.2022
 *      Author: wohlrapp
 */

#include "./SphereInputReader.h"
#include "../utils/ArrayUtils.h"
#include "../utils/ParticleGenerator.h"

#include <string>
#include <fstream>
#include <sstream>

void SphereInputReader::readInput(ProgramParameters &programParameters, const char *filename)
{
    // Variables to read in
    std::array<double, 3> center;
    int r;
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

        for (auto &xi : center)
        {
            datastream >> xi;
        }
        datastream.clear();

        // get next line which contains the n
        getline(input_file, tmp_string);
        getLogicLogger()->info("Read line: {}", tmp_string);
        datastream.str(tmp_string);

        datastream >> r;
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
    std::unique_ptr<Sphere> sphere = std::make_unique<Sphere>(Sphere(center, r, h, m, v, epsilon, sigma, type, fixed_bool));
    ParticleGenerator::generateSphere(*particleContainer, *sphere, programParameters.getDimension());
}