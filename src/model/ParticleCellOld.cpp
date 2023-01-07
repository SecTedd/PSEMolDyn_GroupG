/*
 * ParticleCell.cpp
 *
 * Created: 02.12.2022
 * Author:  marquardt
 */

#include "ParticleCellOld.h"
#include "../utils/ArrayUtils.h"
#include <iostream>

ParticleCellOld::ParticleCellOld(CellType type, std::array<BoundaryCondition, 6> boundaries)
{
    _type = type;
    _boundaries = boundaries;
    _particles.reset(new std::vector<Particle *>);

    _memoryLogger = spdlog::get("memory_logger");
    _memoryLogger->info("ParticleCell generated!");
    _simulationLogger = spdlog::get("simulation_logger");
}

ParticleCellOld::~ParticleCellOld()
{
    _memoryLogger->info("ParticleCell destructed!");
}

const void ParticleCellOld::insertParticle(Particle *p)
{
    _particles->push_back(p);
}

const void ParticleCellOld::iterateParticlePairs(std::function<void(Particle &, Particle &)> f, double cutoff)
{
    // Since we use a vector we can directly access the particles through indexing
    // Because f_ij = -f_ji, we can save (n^2)/2 iterations by starting the inner loop at i+1
    for (long unsigned int i = 0; i < _particles->size(); i++)
    {
        if (!_particles->at(i)->getInvalid() && !_particles->at(i)->getHalo())
        {
            for (long unsigned int j = i + 1; j < _particles->size(); j++)
            {   if (_particles->at(i)->getX()[0] == _particles->at(j)->getX()[0] && _particles->at(i)->getX()[1] == _particles->at(j)->getX()[1] && _particles->at(i)->getX()[2] == _particles->at(j)->getX()[2]) {
                    _simulationLogger->debug("Particles at same position: " + _particles->at(i)->toString() + _particles->at(j)->toString());
                }
                if (!_particles->at(j)->getInvalid() && !_particles->at(j)->getHalo() && ArrayUtils::L2Norm(_particles->at(i)->getX() - _particles->at(j)->getX()) <= cutoff)
                {
                    f(*_particles->at(i), *_particles->at(j));
                }
            }
        }
    }
}

const void ParticleCellOld::clearCell() { _particles->clear(); }

const void ParticleCellOld::reserveMemory(int meanParticles)
{
    _particles->reserve(_particles->size() + meanParticles);
}

std::vector<Particle *> *ParticleCellOld::getCellParticles()
{
    return _particles.get();
}

const void ParticleCellOld::removeInvalid()
{ 
    _particles->erase(std::remove_if(_particles->begin(), _particles->end(), [](Particle *p)
                                    { return p->getInvalid() || p->getHalo(); }),
                        _particles->end());  
}

const CellType ParticleCellOld::getType() { return _type; }

const std::array<BoundaryCondition, 6> &ParticleCellOld::getBoundaries() { return _boundaries; }

const int ParticleCellOld::size() { return _particles->size(); }

const std::vector<int> &ParticleCellOld::getNeighbours() { return _neighbours; }

const void ParticleCellOld::setNeighbours(std::vector<int> &neighbours) { _neighbours = neighbours; }

const std::vector<int> &ParticleCellOld::getHaloNeighbours() { return _haloNeighbours; }

const void ParticleCellOld::setHaloNeighbours(std::vector<int> &haloNeighbours) { _haloNeighbours = haloNeighbours; }
