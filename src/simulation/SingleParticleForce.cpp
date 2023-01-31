/*
 *  SingleParticleForce.cpp
 *
 *  Created on: 15.12.2022
 *      Author: wohlrapp
 */

#include "SingleParticleForce.h"

SingleParticleForce::SingleParticleForce(std::array<double, 3> force)
{
    this->force = force;
    _memoryLogger = spdlog::get("memory_logger");
    _memoryLogger->info("ForceCalculation generated!");
}

SingleParticleForce::~SingleParticleForce()
{
    _memoryLogger->info("ForcCalculation destructed!");
}