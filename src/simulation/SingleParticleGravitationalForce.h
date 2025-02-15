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
class SingleParticleGravitationalForce : public SingleParticleForce
{
public:
    /**
     * @param force acting on the particles
    */
    SingleParticleGravitationalForce(std::array<double, 3> force); 
    /**
     * @brief calculates the gravitational force acting on a single particle
     *
     * @param particleContainer container for the particles for which the force will be calculated
     * @param time the current time within the simulation
     */
    void calculateForce(ParticleContainer &particleContainer, double time) override;
};