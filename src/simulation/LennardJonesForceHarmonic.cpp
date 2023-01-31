/*
 * LennardJonesSimulation.cpp
 *
 * Created: 16.11.2022
 * Author:  marquardt
 */

#include "LennardJonesForceHarmonic.h"
#include "../model/ParticleContainer.h"
#include "../utils/ArrayUtils.h"

#include <iostream>
#include <math.h>

LennardJonesForceHarmonic::LennardJonesForceHarmonic()
{
    _logicLogger = spdlog::get("simulation_logger");
}

void LennardJonesForceHarmonic::calculateForce(ParticleContainer &particleContainer)
{
    // first we iterate over each particle once to initialize new force vector to zero
    std::function<void(Particle &)> forceInitializationIteration = [](Particle &p1)
    {
        p1.setOldF(p1.getF());
        p1.setF({0.0, 0.0, 0.0});
    };

    particleContainer.iterateParticles(forceInitializationIteration, false);

    std::function<void(Particle &)> harmonicPotential = [&](Particle &p1)
    {
        double stiffness = p1.getStiffness();
        double averageBondLength = p1.getAverageBondLength();

        auto &particles = particleContainer.getActiveParticles();

        for (int index : p1.getParallelNeighbours())
        {
            double distance = ArrayUtils::L2Norm(p1.getX() - particles[index].getX());
            auto f_ij = (stiffness * (distance - averageBondLength) / distance) * (particles[index].getX() - p1.getX());
            auto f_ji = -1 * f_ij;
            p1.addF(f_ij);
            particles[index].addF(f_ji);
        }

        auto sqrt = std::sqrt(2);

        for (int index : p1.getDiagonalNeighbours())
        {
            double distance = ArrayUtils::L2Norm(p1.getX() - particles[index].getX());
            auto f_ij = (stiffness * (distance - sqrt * averageBondLength) / distance) * (particles[index].getX() - p1.getX());
            auto f_ji = -1 * f_ij;
            p1.addF(f_ij);
            particles[index].addF(f_ji);
        }
    };
    particleContainer.iterateParticles(harmonicPotential, false);

    // in the second step we calculate the forces between pairs of particles according to the Lennard-Jones formula
    std::function<void(Particle &, Particle &)> forceCalculationIteration = [this](Particle &p1, Particle &p2)
    {
        double sigma = (p1.getSigma() + p2.getSigma()) / 2;
        double distance = ArrayUtils::L2Norm(p1.getX() - p2.getX());

        if (distance < std::pow(2, 1.0 / 6) * sigma)
        {
            double epsilon;
            if (p1.getEpsilon() == p2.getEpsilon())
                epsilon = p1.getEpsilon();
            else
                epsilon = sqrt(p1.getEpsilon() * p2.getEpsilon());

            // Reduce number of operation by reusing previous results
            double pow1 = sigma / distance;
            double pow2 = pow1 * pow1;
            double pow4 = pow2 * pow2;
            double pow6 = pow4 * pow2;
            double pow12 = pow6 * pow6;

            // Lennard-Jones force
            std::array<double, 3> f_ij = (-24 * epsilon / pow(distance, 2)) * (pow6 - 2 * pow12) * (p1.getX() - p2.getX());

            std::array<double, 3> f_ji = -1 * f_ij;

            p1.addF(f_ij);
            p2.addF(f_ji);
        }
    };

    particleContainer.iterateParticleInteractions(forceCalculationIteration);
}