/*
 *  SingleParticleGraviationalForce.h
 *
 *  Created on: 15.12.2022
 *      Author: wohlrapp
 */

#pragma once

#include "SingleParticleForce.h"

/**
 * @brief Class calculating the gravitational force acting on one particle
 */
class TemporalSingleParticleForce : public SingleParticleForce
{
private: 
    double endTime;
    std::vector<int> particles;
public:
    /**
     * @param force acting on the particles 
     * @param endTime the time until when the force acts on the particles 
     * @param particles the indices of the affected particles
    */
    TemporalSingleParticleForce(std::array<double, 3> force, double endTime, std::vector<int> particles); 
    /**
     * @brief calculates the gravitational force acting on a single particle
     *
     * @param particleContainer container for the particles for which the force will be calculated
     * @param g_grav the gravitational force acting on a single particle
     */
    void calculateForce(ParticleContainer &particleContainer, double time) override;
};