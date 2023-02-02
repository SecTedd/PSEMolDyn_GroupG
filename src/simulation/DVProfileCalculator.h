/*
 *  DVProfileCalculator.h
 *
 *  Created on: 25.01.2023
 *      Author: borisov
 */

#pragma once

#include "../model/Particle.h"
#include "../model/ParticleContainer.h"
#include "spdlog/spdlog.h"

/**
 * @brief class that calculates the density and velocity profile along x axis
 */
class DVProfileCalculator {
private:
    // subdivision of the x axis (num bins)
    int n;

    // averages per bin
    std::vector<int> avg;

    // domain of the simulation
    std::array<double, 3> domain;

    // particles of the observed system
    std::shared_ptr<ParticleContainer> particleContainer;
    
    // a speedlog logger which logs the logic flow of the simulation
    std::shared_ptr<spdlog::logger> _logicLogger;
    
    //a speedlog logger which logs construction and destruction of particles
    std::shared_ptr<spdlog::logger> _memoryLogger;

public:

    /**
     * @brief Construct a new DVProfileCalculator object
     * 
     * @param particleContainer particles of the observed system
     * @param n number of bins the x axis will be divided into
     * @param domain domain of simulation in each dimension
     */
    DVProfileCalculator(std::shared_ptr<ParticleContainer> particleContainer, int n, std::array<double, 3> domain);

    // Destructor
    ~DVProfileCalculator();

    /**
     * @brief calculates the velocity and density profile along x-axis
     * 
     * @return std::vector<int> the number of particles per bin
     */
    const std::vector<int> calculate();

    // Getters
    const int getN();
    const std::vector<int> getAvg();

    // Setters
    const void setN(int n);
    const void setAvg(std::vector<int> avg);
};
