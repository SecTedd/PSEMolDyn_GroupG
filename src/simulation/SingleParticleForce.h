/*
 *  SingleParticleForce.h
 *
 *  Created on: 15.12.2022
 *      Author: wohlrapp
 */

#pragma once

#include "../model/ParticleContainer.h"

/**
 * @brief abstract class for force calculation acting on single particles
 */
class SingleParticleForce
{
private:
    /**
     * a speedlog logger which logs construction and destruction of particles
     */
    std::shared_ptr<spdlog::logger> _memoryLogger;

    std::array<double, 3> force;

public:
    /**
     * @brief calculates the force acting on every single particle
     *
     * @param force the force that acts on each particle
     */
    SingleParticleForce(std::array<double, 3> force);

    virtual ~SingleParticleForce() = 0;

    /**
     * @brief calculates the force acting on a single particle
     *
     * @param particleContainer container for the particles for which the force will be calculated
     * @param time the current time when the function is called
     */
    virtual void calculateForce(ParticleContainer &particleContainer, double time) = 0;

    const std::array<double, 3> getForce() const { return force; }
};