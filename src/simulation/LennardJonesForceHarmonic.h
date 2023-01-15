/*
 * LennardJonesSimulation.h
 *
 * Created: 16.11.2022
 * Author:  marquardt
 */

#pragma once

#include "InterParticleForce.h"

/**
 * @brief Provides the force calculation according to Lennard-Jones potential
 */
class LennardJonesForceHarmonic : public InterParticleForce
{
private:
    std::shared_ptr<spdlog::logger> _logicLogger;

    double stiffness; /// stiffness of membrane

    double averageBondLength; /// average bond length of membrane

public:
    /**
     * @brief Calculate the force for all particles
     *
     * @param particleContainer The container which contains the particles
     */
    void calculateForce(ParticleContainer &particleContainer) override;

    LennardJonesForceHarmonic(double stiffness, double averageBondLength);
};