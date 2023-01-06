
#include "LinkedCellParticleContainer.h"
#include "../utils/PContainer.h"
#include <cmath>

LinkedCellParticleContainer::LinkedCellParticleContainer(double cutoff, std::array<double, 3> &domain, std::array<BoundaryCondition, 6> &domainBoundaries)
{
    _memoryLogger = spdlog::get("memory_logger");
    _memoryLogger->info("LinkedCellParticleContainer generated!");
    _simulationLogger = spdlog::get("simulation_logger");

    domain = domain;
    cutoff = cutoff;

    initializeCells(domainBoundaries);
}

LinkedCellParticleContainer::~LinkedCellParticleContainer()
{
    _memoryLogger->info("LinkedCellParticleContainer destructed!");
}

const void LinkedCellParticleContainer::initializeCells(std::array<BoundaryCondition, 6> &domainBoundaries)
{
    // smallest cellsize bigger than cutoff plus boundaries
    int numberOfXCells = static_cast<int>(std::floor(domain[0] / _cutoff));
    int numberOfYCells = static_cast<int>(std::floor(domain[1] / _cutoff));
    int numberOfZCells = static_cast<int>(std::floor(domain[2] / _cutoff));

    // for 2D we need exactly one layer, which is not guaranteed above
    if (numberOfZCells == 0)
        numberOfZCells = 1;

    // add surrounding halo cells
    numCells = {numberOfXCells + 2, numberOfYCells + 2, numberOfZCells + 2};

    double sizeX = domain[0] / numberOfXCells;
    double sizeY = domain[1] / numberOfYCells;
    double sizeZ = domain[2] / numberOfZCells;

    // cell size has to be at least 1, otherwise we get problems with division by 0
    if (sizeZ == 0)
    {
        sizeZ = 1;
    }

    cellSize = {sizeX, sizeY, sizeZ};

    // now we need to initialize the cells
    cells.reserve(numCells[0] * numCells[1] * numCells[2]);

    for (int i = 0; i < numCells[0] * numCells[1] * numCells[2]; i++)
    {
        CellType ct = CellType::InnerCell;
        BoundaryCondition of = BoundaryCondition::Outflow;

        std::array<BoundaryCondition, 6> boundaries = {of, of, of, of, of, of};
        std::array<int, 3> index = PContainer::convert1DTo3D(i, numCells);

        // halo cell
        if (index[0] == 0 && domainBoundaries[0] == BoundaryCondition::Periodic || index[0] == numCells[0] - 1 && domainBoundaries[1] == BoundaryCondition::Periodic || index[1] == 0 && domainBoundaries[2] == BoundaryCondition::Periodic || index[1] == numCells[1] - 1 && domainBoundaries[3] == BoundaryCondition::Periodic || index[2] == 0 && domainBoundaries[4] == BoundaryCondition::Periodic || index[2] == numCells[2] - 1 && domainBoundaries[5] == BoundaryCondition::Periodic)
        {
            ct = CellType::PeriodicHaloCell;
        }
        else if (index[0] == 0 || index[0] == numCells[0] - 1 || index[1] == 0 || index[1] == numCells[1] - 1 || index[2] == 0 || index[2] == numCells[2] - 1)
        {
            ct = CellType::HaloCell;
        }

        // cell at left boundary
        if (index[0] == 1)
        {
            boundaries[0] = domainBoundaries[0];
        }
        // cell at right boundary
        if (index[0] == numCells[0] - 2)
        {
            boundaries[1] = domainBoundaries[1];
        }
        // cell at bottom boundary
        if (index[1] == 1)
        {
            boundaries[2] = domainBoundaries[2];
        }
        // cell at top boundary
        if (index[1] == numCells[1] - 2)
        {
            boundaries[3] = domainBoundaries[3];
        }
        // cell at front boundary
        if (index[2] == 1)
        {
            boundaries[4] = domainBoundaries[4];
        }
        // cell at back boundary
        if (index[2] == numCells[2] - 2)
        {
            boundaries[5] = domainBoundaries[5];
        }

        cells.emplace_back(ct, boundaries);
    }
}

const int LinkedCellParticleContainer::computeCellIdx(Particle &p)
{
    //+1 because of halo cells
    int cell_idx_x = static_cast<int>(std::floor(p.getX()[0] / cellSize[0])) + 1;
    int cell_idx_y = static_cast<int>(std::floor(p.getX()[1] / cellSize[1])) + 1;
    int cell_idx_z = static_cast<int>(std::floor(p.getX()[2] / cellSize[2])) + 1;

    return cell_idx_z * numCells[1] * numCells[0] + cell_idx_y * numCells[0] + cell_idx_x;
}

std::array<double,3> LinkedCellParticleContainer::mirroredPosition(std::array<double, 3> position){
    return {fmod(position[0]+domain[0], domain[0]), fmod(position[1]+domain[1], domain[1]), fmod(position[2]+domain[2], domain[2])};
}






const void LinkedCellParticleContainer::addParticle(std::array<double, 3> &x, std::array<double, 3> &v, double &m, double &epsilon, double &sigma)
{
    if (x[0] >= 0 && x[0] < domain[0] && x[1] >= 0 && x[1] < domain[1] && x[2] >= 0 && x[2] < domain[2])
    {

        activeParticles.emplace_back(x, v, m, epsilon, sigma);
        rebuildCells();
    }
}

const void LinkedCellParticleContainer::addParticle(std::array<double, 3> &x, std::array<double, 3> &v, double &m, double &epsilon, double &sigma, int &type)
{
    if (x[0] >= 0 && x[0] < domain[0] && x[1] >= 0 && x[1] < domain[1] && x[2] >= 0 && x[2] < domain[2])
    {

        activeParticles.emplace_back(x, v, m, epsilon, sigma);
        rebuildCells();
    }
}

const void LinkedCellParticleContainer::addParticle(std::array<double, 3> &x, std::array<double, 3> &v, std::array<double, 3> &f, std::array<double, 3> &old_f, double &m, double &epsilon, double &sigma, int &type)
{
    if (x[0] >= 0 && x[0] < domain[0] && x[1] >= 0 && x[1] < domain[1] && x[2] >= 0 && x[2] < domain[2])
    {

        activeParticles.emplace_back(x, v, m, epsilon, sigma);
        rebuildCells();
    }
}