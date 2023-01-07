
#include "LinkedCellParticleContainer.h"
#include "../utils/PContainer.h"

#include <iostream>
#include <cmath>

LinkedCellParticleContainer::LinkedCellParticleContainer(double cutoff, std::array<double, 3> &domain, std::array<BoundaryCondition, 6> &domainBoundaries)
{
    _memoryLogger = spdlog::get("memory_logger");
    _memoryLogger->info("LinkedCellParticleContainer generated!");
    _simulationLogger = spdlog::get("simulation_logger");

    this->domain = domain;
    this->cutoff = cutoff;

    initializeCells(domainBoundaries);
}

LinkedCellParticleContainer::~LinkedCellParticleContainer()
{
    _memoryLogger->info("LinkedCellParticleContainer destructed!");
}

const void LinkedCellParticleContainer::initializeCells(std::array<BoundaryCondition, 6> &domainBoundaries)
{
    // smallest cellsize bigger than cutoff plus boundaries
    int numberOfXCells = static_cast<int>(std::floor(domain[0] / cutoff));
    int numberOfYCells = static_cast<int>(std::floor(domain[1] / cutoff));
    int numberOfZCells = static_cast<int>(std::floor(domain[2] / cutoff));

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
        sizeZ = 1;

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
        if ((index[0] == 0 && domainBoundaries[0] == BoundaryCondition::Periodic) || (index[0] == numCells[0] - 1 && domainBoundaries[1] == BoundaryCondition::Periodic) || (index[1] == 0 && domainBoundaries[2] == BoundaryCondition::Periodic) || (index[1] == numCells[1] - 1 && domainBoundaries[3] == BoundaryCondition::Periodic) || (index[2] == 0 && domainBoundaries[4] == BoundaryCondition::Periodic) || (index[2] == numCells[2] - 1 && domainBoundaries[5] == BoundaryCondition::Periodic))
        {
            ct = CellType::PeriodicHaloCell;
        }
        else if (index[0] == 0 || index[0] == numCells[0] - 1 || index[1] == 0 || index[1] == numCells[1] - 1 || index[2] == 0 || index[2] == numCells[2] - 1)
        {
            ct = CellType::HaloCell;
        }
        else
        {
            // cell at left boundary
            if (index[0] == 1)
            {
                boundaries[0] = domainBoundaries[0];
                ct = CellType::BoundaryCell;
            }
            // cell at right boundary
            if (index[0] == numCells[0] - 2)
            {
                boundaries[1] = domainBoundaries[1];
                ct = CellType::BoundaryCell;
            }
            // cell at bottom boundary
            if (index[1] == 1)
            {
                boundaries[2] = domainBoundaries[2];
                ct = CellType::BoundaryCell;
            }
            // cell at top boundary
            if (index[1] == numCells[1] - 2)
            {
                boundaries[3] = domainBoundaries[3];
                ct = CellType::BoundaryCell;
            }
            // cell at front boundary
            if (index[2] == 1)
            {
                boundaries[4] = domainBoundaries[4];
                ct = CellType::BoundaryCell;
            }
            // cell at back boundary
            if (index[2] == numCells[2] - 2)
            {
                boundaries[5] = domainBoundaries[5];
                ct = CellType::BoundaryCell;
            }
        }
        cells.emplace_back(ct, boundaries);
    }

    for (long unsigned int i = 0; i < cells.size(); i++)
    {
        auto neighbours = PContainer::getDomainNeighboursNewton(i, numCells);
        cells[i].setNeighbours(neighbours);
    }
}

const int LinkedCellParticleContainer::computeCellIdx(Particle &p)
{
    //+1 because of halo cells
    int cell_idx_x = static_cast<int>(std::floor(p.getX()[0] / cellSize[0])) + 1;
    int cell_idx_y = static_cast<int>(std::floor(p.getX()[1] / cellSize[1])) + 1;
    int cell_idx_z = static_cast<int>(std::floor(p.getX()[2] / cellSize[2])) + 1;

    int cellIdx = cell_idx_z * numCells[1] * numCells[0] + cell_idx_y * numCells[0] + cell_idx_x;

    return cellIdx;
}

const void LinkedCellParticleContainer::updateCells()
{
    for (auto &cell : cells)
    {
        if (cell.getType() == CellType::BoundaryCell || cell.getType() == CellType::InnerCell)
            cell.removeInvalid(&activeParticles);
    }

    for (long unsigned int i = 0; i < activeParticles.size(); i++)
    {
        if (activeParticles[i].getInvalid())
        {
            int cellIndex = computeCellIdx(activeParticles[i]);
            cells[cellIndex].insertParticleIndex(i);
            activeParticles[i].setInvalid(false);
        }
    }
}

const void LinkedCellParticleContainer::rebuildCells()
{
    for (auto &cell : cells)
    {
        cell.clearCell();
    }

    // first we need to remove all halo particles
    activeParticles.erase(std::remove_if(activeParticles.begin(), activeParticles.end(), [](Particle p)
                                         { return p.getHalo(); }),
                          activeParticles.end());

    // then rebuild the cells
    for (long unsigned int i = 0; i < activeParticles.size(); i++)
    {
        int cellIndex = computeCellIdx(activeParticles[i]);
        cells[cellIndex].insertParticleIndex(i);
        activeParticles[i].setInvalid(false);
    }
}

std::array<double, 3> LinkedCellParticleContainer::mirroredPosition(std::array<double, 3> position)
{
    return {fmod(position[0] + domain[0], domain[0]), fmod(position[1] + domain[1], domain[1]), fmod(position[2] + domain[2], domain[2])};
}

const void LinkedCellParticleContainer::iterateParticleInteractions(std::function<void(Particle &, Particle &)> f)
{
    // first we need to clear the halo
    clearHalo();

    // then add particles to halo or add the reflecting force
    for (long unsigned int i = 0; i < cells.size(); i++)
    {
        if (cells[i].getType() == CellType::BoundaryCell)
        {
            // handles reflecting boundaries
            reflectingBoundary(i, f);
            // handles periodic boundaries
            periodicBoundary(i);
        }
    }

    // finally calculate interactions
    for (auto &cell : cells)
    {
        // inner forces don't need to be calculated for halo cells
        if (cell.getType() == CellType::InnerCell || cell.getType() == CellType::BoundaryCell)
        {
            std::vector<int> *particleIndices = cell.getCellParticleIndices();
            // iterate within the cell first
            for (long unsigned int i = 0; i < particleIndices->size(); i++)
            {
                for (long unsigned int j = i + 1; j < particleIndices->size(); j++)
                {
                    if (ArrayUtils::L2Norm(activeParticles[i].getX() - activeParticles[j].getX()) <= cutoff)
                    {
                        f(activeParticles[i], activeParticles[j]);
                    }
                }
            }
        }

        // then iterate interactions between particles from different cells
        for (auto particleIndex : *cell.getCellParticleIndices())
        {
            for (auto &neighbouringCellIndex : cell.getNeighbours())
            {
                for (auto interactingParticleIndex : *cells[neighbouringCellIndex].getCellParticleIndices())
                {
                    // figure out if cells are halo or domain cells
                    Particle *currentParticle;
                    if (cell.getType() == CellType::HaloCell || cell.getType() == CellType::PeriodicHaloCell)
                    {
                        currentParticle = &haloParticles[particleIndex];
                    }
                    else
                    {
                        currentParticle = &activeParticles[particleIndex];
                    }

                    Particle *interactingParticle;
                    if (cells[neighbouringCellIndex].getType() == CellType::HaloCell || cells[neighbouringCellIndex].getType() == CellType::PeriodicHaloCell)
                    {
                        interactingParticle = &haloParticles[interactingParticleIndex];
                    }
                    else
                    {
                        interactingParticle = &activeParticles[interactingParticleIndex];
                    }
                    if (ArrayUtils::L2Norm(currentParticle->getX() - interactingParticle->getX()) <= cutoff)
                    {
                        f(*currentParticle, *interactingParticle);
                    }
                }
            }
        }
    }
}

const void LinkedCellParticleContainer::iterateParticles(std::function<void(Particle &)> f, bool calcX)
{
    // no deletions in the active particle vector, meaning it suffices to delete invalid particles from the cells and add the particles elsewhere
    bool cellUpdate = false;
    // deletions in the active particle vector, meaning indices within the cells are rendered invalid
    bool cellRebuild = false;

    for (auto &particle : activeParticles)
    {

        f(particle);

        // check if particle moved out of cell
        int cellIndex = computeCellIdx(particle);

        if (cellIndex != particle.getCellIdx())
        {
            // particle can either be in a regular halo cell which means outflow, a periodic halo which means it needs to be mirrored or in an inner cell which means the index has to be changed
            if (cells[cellIndex].getType() == CellType::InnerCell || cells[cellIndex].getType() == CellType::BoundaryCell)
            {
                // only change in cell index and then a rebuild required
                particle.setCellIdx(cellIndex);
                particle.setInvalid(true);
                cellUpdate = true;
            }
            else if (cells[cellIndex].getType() == CellType::HaloCell)
            {
                particle.setCellIdx(cellIndex);
                particle.setHalo(true);
                cellRebuild = true;
            }
            else if (cells[cellIndex].getType() == CellType::PeriodicHaloCell)
            {
                particle.setX(mirroredPosition(particle.getX()));
                particle.setInvalid(true);
                cellUpdate = true;
            }
        }
    }

    // reorganising cell structure
    if (cellRebuild)
        rebuildCells();
    else if (cellUpdate)
        updateCells();
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

const void LinkedCellParticleContainer::reflectingBoundary(int cellIndex, std::function<void(Particle &, Particle &)> f)
{
    auto boundaries = cells[cellIndex].getBoundaries();

    Particle counterParticle = Particle();

    for (long unsigned int b = 0; b < boundaries.size(); b++)
    {
        if (boundaries[b] == BoundaryCondition::Reflecting)
        {
            std::cout << "reflecting boundary at: " << b << std::endl;
            for (auto particleIndex : *cells[cellIndex].getCellParticleIndices())
            {
                double reflectingDistance = std::pow(2, 1.0 / 6) * activeParticles[particleIndex].getSigma() / 2;
                counterParticle.setEpsilon(activeParticles[particleIndex].getEpsilon());
                counterParticle.setSigma(activeParticles[particleIndex].getSigma());
                switch (b)
                {
                // left boundary, x-coordinate is 0
                case 0:
                    if (activeParticles[particleIndex].getX()[0] <= reflectingDistance)
                    {
                        std::array<double, 3> counterX = {(-1) * activeParticles[particleIndex].getX()[0], activeParticles[particleIndex].getX()[1], activeParticles[particleIndex].getX()[2]};
                        counterParticle.setX(counterX);
                        f(activeParticles[particleIndex], counterParticle);
                    }
                    break;
                // right boundary, x-coordinate is domain size in x-direction
                case 1:
                    if (domain[0] - activeParticles[particleIndex].getX()[0] <= reflectingDistance)
                    {
                        std::array<double, 3> counterX = {domain[0] + domain[0] - activeParticles[particleIndex].getX()[0], activeParticles[particleIndex].getX()[1], activeParticles[particleIndex].getX()[2]};
                        counterParticle.setX(counterX);
                        f(activeParticles[particleIndex], counterParticle);
                    }
                    break;
                // lower boundary, y-coordinate is 0
                case 2:
                    if (activeParticles[particleIndex].getX()[1] <= reflectingDistance)
                    {
                        std::array<double, 3> counterX = {activeParticles[particleIndex].getX()[0], (-1) * activeParticles[particleIndex].getX()[1], activeParticles[particleIndex].getX()[2]};
                        counterParticle.setX(counterX);
                        f(activeParticles[particleIndex], counterParticle);
                    }
                    break;
                // upper boundary, y-coordinate is domain size in y-direction
                case 3:
                    if (domain[1] - activeParticles[particleIndex].getX()[1] <= reflectingDistance)
                    {
                        std::array<double, 3> counterX = {activeParticles[particleIndex].getX()[0], domain[1] + domain[1] - activeParticles[particleIndex].getX()[1], activeParticles[particleIndex].getX()[2]};
                        counterParticle.setX(counterX);
                        f(activeParticles[particleIndex], counterParticle);
                    }
                    break;
                // front boundary, z-coordinate is 0
                case 4:
                    if (activeParticles[particleIndex].getX()[2] <= reflectingDistance)
                    {
                        std::array<double, 3> counterX = {activeParticles[particleIndex].getX()[0], activeParticles[particleIndex].getX()[1], (-1) * activeParticles[particleIndex].getX()[2]};
                        counterParticle.setX(counterX);
                        f(activeParticles[particleIndex], counterParticle);
                    }
                    break;
                // back boundary, z-coordinate is domain size in z-direction
                case 5:
                    if (domain[2] - activeParticles[particleIndex].getX()[2] <= reflectingDistance)
                    {
                        std::array<double, 3> counterX = {activeParticles[particleIndex].getX()[0], activeParticles[particleIndex].getX()[1], domain[2] + domain[2] - activeParticles[particleIndex].getX()[2]};
                        counterParticle.setX(counterX);
                        f(activeParticles[particleIndex], counterParticle);
                    }
                    break;
                }
            }
        }
    }
}

const void LinkedCellParticleContainer::periodicBoundary(int cellIndex)
{
    // create a list with every boundary that is periodic
    std::vector<int> periodicBoundaries;

    for (long unsigned int b = 0; b < cells[cellIndex].getBoundaries().size(); b++)
    {
        if (cells[cellIndex].getBoundaries()[b] == BoundaryCondition::Periodic)
            periodicBoundaries.push_back(b);
    }
    if (periodicBoundaries.size() > 0)
    {
        // function which returns indices where cell is supposed to be mirrored to
        auto mirroringOffsets = PContainer::getMirroringOffsets(cellIndex, periodicBoundaries, numCells, cellSize);

        // mirror current cell to every relevant cell
        for (auto mirroringOffset : mirroringOffsets)
        {
            for (auto particleIndex : *cells[cellIndex].getCellParticleIndices())
            {
                // std::array<double, 3> &x, std::array<double, 3> &v, double &m, double &epsilon, double &sigma
                std::array<double, 3> newX = {activeParticles[particleIndex].getX()[0] + mirroringOffset[0], activeParticles[particleIndex].getX()[1] + mirroringOffset[1], activeParticles[particleIndex].getX()[2] + mirroringOffset[2]};
                haloParticles.emplace_back(newX, activeParticles[particleIndex].getV(), activeParticles[particleIndex].getM(), activeParticles[particleIndex].getEpsilon(), activeParticles[particleIndex].getSigma());
            }
        }
    }
}

const void LinkedCellParticleContainer::clearHalo()
{
    haloParticles.clear();

    for (auto &cell : cells)
    {
        if (cell.getType() == CellType::PeriodicHaloCell || cell.getType() == CellType::HaloCell)
            cell.clearCell();
    }
}

const int LinkedCellParticleContainer::size() const
{
    return activeParticles.size();
}

const void LinkedCellParticleContainer::resetParticles()
{
    activeParticles.clear();
    haloParticles.clear();

    for (auto &cell : cells)
        cell.clearCell();
}

const void LinkedCellParticleContainer::reserveMemoryForParticles(int numberOfParticles)
{
    // reserving extra space before new particles are added. Push_back can then be executed without resizing
    int newLength = numberOfParticles + activeParticles.size();
    activeParticles.reserve(newLength);
    haloParticles.reserve(newLength);
    // reserving extra space in each cell with the mean value of particles per cell
    for (auto cell : cells)
    {
        // for number of mean particles per cell consider only domain cells
        cell.reserveMemory(numberOfParticles / ((numCells[0] - 2) * (numCells[1] - 2) * (numCells[2] - 2)));
    }
    rebuildCells();
}

std::vector<Particle> *LinkedCellParticleContainer::getBoundaryParticles()
{

    auto boundaryParticles = std::make_shared<std::vector<Particle>>();

    for (auto &cell : cells)
    {
        if (cell.getType() == CellType::BoundaryCell)
        {
            for (auto particleIndex : *cell.getCellParticleIndices())
                boundaryParticles->push_back(activeParticles[particleIndex]);
        }
    }

    return boundaryParticles.get();
}

std::vector<ParticleCell> &LinkedCellParticleContainer::getCells() { return cells; }

std::vector<Particle> &LinkedCellParticleContainer::getHaloParticles() { return haloParticles; }

std::vector<Particle> &LinkedCellParticleContainer::getActiveParticles() { return activeParticles; }