/*
 * ParticleCell.h
 *
 * Created: 02.12.2022
 * Author:  marquardt
 */

#pragma once

#include "Particle.h"

#include <vector>

/**
 * @brief Enum for different types of cells
 */
enum class CellType
{
    InnerCell,
    BoundaryCell,
    HaloCell, 
    PeriodicHaloCell
};

/**
 * @brief Enum for different types of boundary conditions
 */
enum class BoundaryCondition
{
    Outflow,
    Reflecting,
    Periodic
};

/**
 * @brief class for cells used in the linked cell algorithm
 */
class ParticleCell
{
private:
    std::shared_ptr<std::vector<int>> particleIndices; //reference to vector of pointers to particles currently in this cell

    std::vector<int> domainNeighbours; // structure to store index of neighboring domain cells with a higher index

    std::vector<int> haloNeighbours; //structure to store all neighbouring halo cells

    CellType type; // type of cell (inner or boundary)

    /**
     * array to store boundary conditions for each cell
     * first two indices for x-direction (left, right)
     * middle two indices for y-direction (bottom, top)
     * last two indices for z-direction (front, back)
     */
    std::array<BoundaryCondition, 6> boundaries;

    std::shared_ptr<spdlog::logger> memoryLogger; // a speedlog logger which logs construction and destruction of particles
    std::shared_ptr<spdlog::logger> simulationLogger;

public:
    /**
     * @brief constructor for the particle cell, initializes loggers and variables
     * @param type type of cell (inner, boundary)
     * @param boundaries boundary condition on each border of the cell
     */
    ParticleCell(CellType type, std::array<BoundaryCondition, 6> boundaries);
    ~ParticleCell();

    /**
     * @brief inserts pointer to particle at the end of particle vector
     * @param p pointer to new particle
     */
    const void insertParticleIndex(int index);

    /**
     * @brief removes all particles from cell
     */
    const void clearCell();

    /**
     * @brief reserves extra vector space for given number of particles
     * @param meanParticles mean number of particles per cell, estimate for actual number of particles
     */
    const void reserveMemory(int meanParticles);

    /**
     * @brief removes particle pointers of invalid particles
    */
    const void removeInvalid(std::vector<Particle> *particles);

    /**
     * @brief returns particles of this cell
     * @return pointer to particle vector of this cell
     */
    std::vector<int> *getCellParticleIndices();

    const std::array<BoundaryCondition, 6> &getBoundaries();

    const CellType getType();

    const int size();

    const std::string toString();

    const std::vector<int> &getDomainNeighbours();

    const std::vector<int> &getHaloNeighbours();

    void setDomainNeighbours(std::vector<int> &neighbours);

    void setHaloNeighbours(std::vector<int> &neighbours);
};