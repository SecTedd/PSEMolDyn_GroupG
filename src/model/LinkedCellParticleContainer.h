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

    std::vector<std::vector<int>> superCells; // clusters together a batch of 4 (2D) or 8 (3D) neighbouring cells

    std::array<double, 3> domain; // domain size in each dimension

    double cutoff; // max. distance of particles where force calculation is applied

    std::array<int, 3> numCells; // number of cells in each dimension

    std::array<double, 3> cellSize; // cell size in each dimension

    std::array<int, 3> numInteractingCells; // number of cells in each dimension relevant for grouping

    int parallel; // 0 for no parallelization, 1 for fork join with cell grouping, 2 for fork join with supercell grouping, 3 for tasks

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
     * @brief reserves memory for cell groups according to parallelization strategy
     */
    void reserveGroups();

    /**
     * @brief reserves memory for vector of supercells, only called if parallel = 2
     */
    void reserveSuperCells();

    /**
     * @brief fills group & supercell vectors according to parallelization strategy
     */
    void initializeParallelGroups();

    /**
     * @brief computes group given supercell belongs to
     * @param cellIdx index of supercell
     * @return group index of supercell
     */
    const int computeSupercellGroup(int cellIdx);

    /**
     * @brief computes group or supercell (if parallel = 2) the given cell belongs to according to parallelization strategy
     * @param cellIdx index of cell
     * @return group or supercell index of cell
     */
    const int computeCellGroup(int cellIdx);

    /**
     * @brief computes particle interactions within one cell implementing Newton's 3rd law
     * @param i index of cell
     * @param f particle interaction function
     */
    inline void intraCellInteraction(int i, std::function<void(Particle &, Particle &)> f);

    /**
     * @brief computes particle interactions between particles of two cells
     * @param i index of first cell
     * @param j index of second cell
     * @param f particle interaction function
     */
    inline void interCellInteraction(int i, int j, std::function<void(Particle &, Particle &)> f);

    /**
     * @brief parallelizes particle interactions according to fork join model
     * @param f particle interaction function
     */
    void forkJoin(std::function<void(Particle &, Particle &)> f);

    /**
     * @brief parallelizes neighbouring cell interactions through sequential groups
     * @param f particle interaction function
     */
    void directCellInteraction(std::function<void(Particle &, Particle &)> f);

    /**
     * @brief parallelizes neighbouring cell interactions through supercells and sequential groups
     * @param f particle interaction function
     */
    void nestedCellInteraction(std::function<void(Particle &, Particle &)> f);

    /**
     * @brief parallelizes particle interactions according to task model
     * @param f particle interaction function
     */
    void taskModel(std::function<void(Particle &, Particle &)> f);

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
     */
    const void initializeCells(std::array<BoundaryCondition, 6> &domainBoundaries);

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

    std::vector<std::vector<int>> &getSuperCells();

    void setParallel(int parallel);

    const int getParallel();
};
