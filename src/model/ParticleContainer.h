/*
 * ParticleContainer.h
 *
 * Created: 3.12.2022
 * Author:  marquardt
 */

#pragma once

#include "Particle.h"

/**
 * @brief abstract class for Particle Container
 */
class ParticleContainer
{

protected:
    /**
     * speedlog logger which logs construction and destruction of particles
     */
    std::shared_ptr<spdlog::logger> _memoryLogger;

    /**
     * speedlog logger which logs information about particle behavior (e.g. leaving cell/domain boundaries)
     */
    std::shared_ptr<spdlog::logger> _simulationLogger;

public:
    virtual ~ParticleContainer(){};

    /**
     * @brief Iterates over all active particles (inside domain) and applies the function f
     * @param f A lambda function applied for every particle
     */
    virtual const void iterateParticles(std::function<void(Particle &)> f, bool calcX) = 0;

    /**
     * @brief Computes interaction between two particles or a particle and a border
     * @param f A lambda function applied between particles or a particle and a border
     */
    virtual const void iterateParticleInteractions(std::function<void(Particle &, Particle &)> f) = 0;

    /**
     * @brief Creates a new particle and adds it to the vector
     * @param x The position array of the particle
     * @param v The velocity array of the particle
     * @param m The mass of the particle
     * @param epsilon The epsilon of the particle
     * @param sigma The sigma of the particle
     */
    virtual const void addParticle(std::array<double, 3> &x, std::array<double, 3> &v, double &m, double &epsilon, double &sigma) = 0;

    /**
     * @brief Creates a new particle and adds it to the vector
     * @param x The position array of the particle
     * @param v The velocity array of the particle
     * @param m The mass of the particle
     * @param type The type of the particle
     * @param epsilon The epsilon of the particle
     * @param sigma The sigma of the particle
     */
    virtual const void addParticle(std::array<double, 3> &x, std::array<double, 3> &v, double &m, double &epsilon, double &sigma, int &type) = 0;

    /**
     * @brief Creates a new particle and adds it to the vector
     * @param x The position array of the particle
     * @param v The velocity array of the particle
     * @param m The mass of the particle
     * @param type The type of the particle
     * @param epsilon The epsilon of the particle
     * @param sigma The sigma of the particle
     * @param fixed The stationary indicator
     */
    virtual const void addParticle(std::array<double, 3> &x, std::array<double, 3> &v, double &m, double &epsilon, double &sigma, bool &fixed, int &type) = 0;


    /**
     * @brief Creates a new particle and adds it to the vector
     * @param x The position array of the particle
     * @param v The velocity array of the particle
     * @param f The force acting on the particle
     * @param old_f The previous force acting on the particle
     * @param m The mass of the particle
     * @param type The type of the particle
     * @param epsilon The epsilon of the particle
     * @param sigma The sigma of the particle
     */
    virtual const void addParticle(std::array<double, 3> &x, std::array<double, 3> &v, std::array<double, 3> &f, std::array<double, 3> &old_f, double &m, double &epsilon, double &sigma, int &type) = 0;

    /**
     * @brief returns the number of active particles
     * @return size of the particle vector
     */
    virtual const int size() const = 0;

    /**
     * @brief deletes all particles (active & halo) from the simulation
     */
    virtual const void resetParticles() = 0;

    /**
     * @brief reserves memory space for given number of particles to avoid constant resizing of the vectors
     * @param numberOfParticles number of additional particles to reserve space for
     */
    virtual const void reserveMemoryForParticles(int numberOfParticles) = 0;

    /**
     * @brief returns vector of active particles (particles inside domain boundaries)
     * @return number of active particles
     */
    virtual std::vector<Particle> &getActiveParticles() = 0;
};