/*
 *  SingleParticleGraviationalForce.cpp
 *
 *  Created on: 15.12.2022
 *      Author: wohlrapp
 */

#include "./SingleParticleGravitationalForce.h"
#include "../model/ParticleContainer.h"
#include "../utils/ArrayUtils.h"

#include <vector>

SingleParticleGravitationalForce::SingleParticleGravitationalForce(std::array<double, 3> force) : SingleParticleForce(force) {}


void SingleParticleGravitationalForce::calculateForce(ParticleContainer &particleContainer, double time)
{
    // we iterate over each particle once and then apply the gravitational force
    std::function<void(Particle &)> forceCalculation = [&](Particle &p1)
    {
        auto force = p1.getM() * getForce(); 
        p1.addF(force);
    };

    particleContainer.iterateParticles(forceCalculation, false);
}