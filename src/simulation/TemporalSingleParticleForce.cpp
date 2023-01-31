/*
 *  SingleParticleGraviationalForce.cpp
 *
 *  Created on: 15.12.2022
 *      Author: wohlrapp
 */

#include "./TemporalSingleParticleForce.h"
#include "../model/ParticleContainer.h"
#include "../utils/ArrayUtils.h"

#include <vector>

TemporalSingleParticleForce::TemporalSingleParticleForce(std::array<double, 3> force, double endTime, std::vector<int> particles) : SingleParticleForce(force), endTime(endTime), particles(particles) {}

void TemporalSingleParticleForce::calculateForce(ParticleContainer &particleContainer, double time)
{
    if (time <= endTime)
    {
        auto &activeParticles = particleContainer.getActiveParticles(); 

        for(auto particle : particles){
            auto force = activeParticles[particle].getM() * getForce();
            activeParticles[particle].addF(force);
        }
    }
}