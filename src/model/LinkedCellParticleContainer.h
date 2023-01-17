/*
 * LinkedCellParticleContainer.h
 *
 * Created: 30.11.2022
 * Author:  marquardt
 */

#pragma once

#include "ParticleCell.h"
#include "ParticleContainer.h"

#include <vector>
#include <array>
#include <functional>

/**
 * @brief Particle Container that implements the linked cell algorithm. The plots below show how the algorithm can improve the performance of the simulation
 */
class LinkedCellParticleContainer : public ParticleContainer
{

private:
    std::vector<Particle> activeParticles; // base vector to store all particles inside the domain

    std::vector<Particle> haloParticles; // base vector to store all halo particles

    std::vector<ParticleCell> cells; // stores all cells

    std::vector<std::vector<int>> cellGroups; // stores indices of grouped cells

    std::array<double, 3> domain; // domain size in each dimension

    double cutoff; // max. distance of particles where force calculation is applied

    std::array<int, 3> numCells; // number of cells in each dimension

    std::array<double, 3> cellSize; // cell size in each dimension

    int parallel; //0 for no parallelization, 1 for first parallel strategy, 2 for second parallel strategy

    /**
     * @brief compute index of cell the given particle belongs to
     * @param p particle
     * @return cell index
     */
    const int computeCellIdx(Particle &p);

    /**
     * @brief called if at least one particle crossed a cell boundary, update particle pointers in cells
     */
    const void updateCells();

    /**
     * @brief called if at least one particle crossed to a halo boundary, leads to restructuring of the cells
     */
    const void rebuildCells();

    /**
     * @brief computes position of ghost particle in halo
     * @param p particle which has to be mirrored
     */
    std::array<double, 3> mirroredPosition(std::array<double, 3> position);

    /**
     * @brief reserves memory for & initializes vector of cell groups according to parallelization strategy
     * @param parallel parallelization strategy which has to be used
    */
    void initializeGroups(int parallel);
    
    /**
     * @brief computes group the given cell belongs to according to parallelization strategy
     * @param cellIdx index of cell
     * @param parallel parallelization strategy 
     * @return group index of cell
    */
    const int computeCellGroup(int cellIdx, int parallel);

    /**
     * @brief computes group the given cell belongs to according to parallelization strategy 1
     * @param cellIdx index of cell
     * @return group index of cell
    */
    const int parallelStrategy1(int cellIdx);

    /**
     * @brief computes group the given cell belongs to according to parallelization strategy 2
     * @param cellIdx index of cell
     * @return group index of cell
    */
    const int parallelStrategy2(int cellIdx);

public:
    LinkedCellParticleContainer(double cutoff, std::array<double, 3> &domain, std::array<BoundaryCondition, 6> &boundaries, int parallel);

    ~LinkedCellParticleContainer() override;

    /**
     * @brief applies given function to particle pairs within the cutoff radius considering each cell and its neighbors, also handles boundary conditions
     * @param f function which is applied to the particle pairs
     */
    const void iterateParticleInteractions(std::function<void(Particle &, Particle &)> f) override;

    /**
     * @brief applies given function to every particle, checks if they cross cell borders
     * @param f function which is applied to the particles
     * @param calcX used to prevent some function calls
     */
    const void iterateParticles(std::function<void(Particle &)> f, bool calcX) override;

    /**
     * @brief adds particle to base vector and its pointer to the cell it belongs to
     * @param x The position array of the particle
     * @param v The velocity array of the particle
     * @param m The mass of the particle
     * @param epsilon The epsilon of the particle
     * @param sigma The sigma of the particle
     */
    const void addParticle(std::array<double, 3> &x, std::array<double, 3> &v, double &m, double &epsilon, double &sigma) override;

    /**
     * @brief Creates a new particle and adds it to the vector
     * @param x The position array of the particle
     * @param v The velocity array of the particle
     * @param m The mass of the particle
     * @param type The type of the particle
     * @param epsilon The epsilon of the particle
     * @param sigma The sigma of the particle
     */
    const void addParticle(std::array<double, 3> &x, std::array<double, 3> &v, double &m, double &epsilon, double &sigma, int &type) override;

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
    const void addParticle(std::array<double, 3> &x, std::array<double, 3> &v, std::array<double, 3> &f, std::array<double, 3> &old_f, double &m, double &epsilon, double &sigma, int &type);

    /**
     * @brief computes number of cells and their size in each dimension, initializes them according to domain boundary conditions
     * @param domainBoundaries the boundaries of the domain
     * @param parallel parallelization strategy which should be used in particleInteractions
     */
    const void initializeCells(std::array<BoundaryCondition, 6> &domainBoundaries, int parallel);

    /**
     * @brief adds reflection force to particles near to the reflecting boundary
     * @param cellIndex the index of the current cell
     * @param f force calculation function which has to be applied at the reflecting boundary
     */
    const void reflectingBoundary(int cellIndex, std::function<void(Particle &, Particle &)> f);

    /**
     * @brief emplaces ghost particles of boundary particles at periodic boundaries in halo
     * @param cellIndex the index of the current cell
     */
    const void periodicBoundary(int cellIndex);

    /**
     * @brief removes halo particles from base vector and clears halo cells, updates particle references in cells
     */
    const void clearHalo();

    /**
     * @brief Returns the number of particles in the simulation
     * @return The size of the particle vector
     */
    const int size() const override;

    /**
     * @brief deletes every particle (active & halo) from simulation, deletes particle pointers in cells
     */
    const void resetParticles() override;

    /**
     * @brief reserves additional memory for particles in basevector and additional memory in cells
     * @param numberOfParticles number of additional particles the space is needed for
     */
    const void reserveMemoryForParticles(int numberOfParticles) override;

    /**
     * @brief a function to get the BoundaryParticles
     * @return the vector of boundary particles
     */
    std::vector<Particle> *getBoundaryParticles();

    std::vector<ParticleCell> &getCells();

    std::vector<Particle> &getHaloParticles();

    std::vector<Particle> &getActiveParticles() override;

    std::vector<std::vector<int>> &getCellGroups();

    void setParallel(int parallel);

    const int getParallel();
};

