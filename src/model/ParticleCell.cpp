/*
 * ParticleCell.cpp
 *
 * Created: 02.12.2022
 * Author:  marquardt
 */

#include "ParticleCell.h"
#include "../utils/ArrayUtils.h"
#include <iostream>

ParticleCell::ParticleCell(CellType type, std::array<BoundaryCondition, 6> boundaries)
{
    this->type = type;
    this->boundaries = boundaries;
    particleIndices.reset(new std::vector<int>);

    memoryLogger = spdlog::get("memory_logger");
    memoryLogger->info("ParticleCell generated!");
    simulationLogger = spdlog::get("simulation_logger");
}

ParticleCell::~ParticleCell()
{
    memoryLogger->info("ParticleCell destructed!");
}

const void ParticleCell::insertParticleIndex(int index)
{
    particleIndices->push_back(index);
}

const void ParticleCell::clearCell() { particleIndices->clear(); }

const void ParticleCell::reserveMemory(int meanParticles)
{
    particleIndices->reserve(particleIndices->size() + meanParticles);
}

std::vector<int> *ParticleCell::getCellParticleIndices()
{
    return particleIndices.get();
}

const void ParticleCell::removeInvalid(std::vector<Particle> *particles)
{
    particleIndices->erase(std::remove_if(particleIndices->begin(), particleIndices->end(), [particles](int particleIndex)
                                     { return particles->at(particleIndex).getInvalid() || particles->at(particleIndex).getHalo(); }),
                      particleIndices->end());
}

const CellType ParticleCell::getType() { return type; }

const std::array<BoundaryCondition, 6> &ParticleCell::getBoundaries() { return boundaries; }

const int ParticleCell::size() { return particleIndices->size(); }

const std::vector<int> &ParticleCell::getNeighbours() { return neighbours; }

void ParticleCell::setNeighbours(std::vector<int> &neighbours) { this->neighbours = neighbours; }
