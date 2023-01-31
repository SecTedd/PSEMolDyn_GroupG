
#include "LinkedCellParticleContainer.h"
#include "../utils/PContainer.h"

#include <iostream>
#include <iomanip>
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
            _memoryLogger->debug("Periodic Halo Cell at " + std::to_string(i));
        }
        else if (index[0] == 0 || index[0] == numCells[0] - 1 || index[1] == 0 || index[1] == numCells[1] - 1 || index[2] == 0 || index[2] == numCells[2] - 1)
        {
            ct = CellType::HaloCell;
            _memoryLogger->debug("Halo Cell at " + std::to_string(i));
        }
        else
        {
            // cell at left boundary
            if (index[0] == 1)
            {
                boundaries[0] = domainBoundaries[0];
                ct = CellType::BoundaryCell;
                _memoryLogger->debug("Left Boundary Cell at " + std::to_string(i));
            }
            // cell at right boundary
            if (index[0] == numCells[0] - 2)
            {
                boundaries[1] = domainBoundaries[1];
                ct = CellType::BoundaryCell;
                _memoryLogger->debug("Right Boundary Cell at " + std::to_string(i));
            }
            // cell at bottom boundary
            if (index[1] == 1)
            {
                boundaries[2] = domainBoundaries[2];
                ct = CellType::BoundaryCell;
                _memoryLogger->debug("Lower Boundary Cell at " + std::to_string(i));
            }
            // cell at top boundary
            if (index[1] == numCells[1] - 2)
            {
                boundaries[3] = domainBoundaries[3];
                ct = CellType::BoundaryCell;
                _memoryLogger->debug("Upper Boundary Cell at " + std::to_string(i));
            }
            // cell at front boundary
            if (index[2] == 1)
            {
                boundaries[4] = domainBoundaries[4];
                ct = CellType::BoundaryCell;
                _memoryLogger->debug("Front Boundary Cell at " + std::to_string(i));
            }
            // cell at back boundary
            if (index[2] == numCells[2] - 2)
            {
                boundaries[5] = domainBoundaries[5];
                ct = CellType::BoundaryCell;
                _memoryLogger->debug("Back Boundary Cell at " + std::to_string(i));
            }
        }
        cells.emplace_back(ct, boundaries);
    }

    for (long unsigned int i = 0; i < cells.size(); i++)
    {
        if (cells[i].getType() == CellType::BoundaryCell || cells[i].getType() == CellType::InnerCell)
        {
            auto domainNeighbours = PContainer::getDomainNeighboursNewton(i, numCells);
            cells[i].setDomainNeighbours(domainNeighbours);
        }
        if (cells[i].getType() == CellType::BoundaryCell)
        {
            std::vector<int> periodicBoundaries;
            for (int b = 0; b < 6; b++)
            {
                if (cells[i].getBoundaries()[b] == BoundaryCondition::Periodic)
                {
                    periodicBoundaries.push_back(b);
                }
            }
            auto haloNeighbours = PContainer::getPeriodicHaloNeighbours(i, numCells, periodicBoundaries);
            cells[i].setPeriodicHaloNeighbours(haloNeighbours);
        }
    }
    _memoryLogger->debug("Cell size: " + std::to_string(sizeX) + ", " + std::to_string(sizeY) + ", " + std::to_string(sizeZ));
    _memoryLogger->debug("Num cells: " + std::to_string(numCells[0]) + ", " + std::to_string(numCells[1]) + ", " + std::to_string(numCells[2]));
}

const int LinkedCellParticleContainer::computeCellIdx(Particle &p)
{
    //+1 because of halo cells
    int cell_idx_x = static_cast<int>(std::floor(p.getX()[0] / cellSize[0])) + 1;
    int cell_idx_y = static_cast<int>(std::floor(p.getX()[1] / cellSize[1])) + 1;
    int cell_idx_z = static_cast<int>(std::floor(p.getX()[2] / cellSize[2])) + 1;
    if (cell_idx_x == numCells[0])
        cell_idx_x = numCells[0] - 1;
    if (cell_idx_y == numCells[1])
        cell_idx_y = numCells[1] - 1;
    if (cell_idx_z == numCells[2])
        cell_idx_z = numCells[2] - 1;
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
        Particle &p = activeParticles[i];
        if (p.getInvalid())
        {
            int cellIndex = computeCellIdx(p);
            cells[cellIndex].insertParticleIndex(i);
            p.setInvalid(false);
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
    activeParticles.erase(std::remove_if(activeParticles.begin(), activeParticles.end(), [](Particle &p)
                                         { return p.getHalo(); }),
                          activeParticles.end());

    // then rebuild the cells
    for (long unsigned int i = 0; i < activeParticles.size(); i++)
    {
        Particle &p = activeParticles[i];
        int cellIndex = computeCellIdx(p);
        cells[cellIndex].insertParticleIndex(i);
        p.setInvalid(false);
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
        // particle interactions don't need to be calculated for halo cells
        if (cell.getType() == CellType::InnerCell || cell.getType() == CellType::BoundaryCell)
        {
            // particle interactions within cell
            std::vector<int> *particleIndices = cell.getCellParticleIndices();
            for (long unsigned int i = 0; i < particleIndices->size(); i++)
            {
                for (long unsigned int j = i + 1; j < particleIndices->size(); j++)
                {
                    Particle &p1 = activeParticles[particleIndices->at(i)];
                    Particle &p2 = activeParticles[particleIndices->at(j)];
                    if (ArrayUtils::L2Norm(p1.getX() - p2.getX()) <= cutoff)
                    {
                        f(p1, p2);
                    }
                }
            }

            std::vector<int> neighbours = cell.getDomainNeighbours();
            // Boundary cells
            if (cell.getType() == CellType::BoundaryCell)
            {
                neighbours.insert(neighbours.end(), cell.getPeriodicHaloNeighbours().begin(), cell.getPeriodicHaloNeighbours().end());
            }

            // particle interactions between different cells
            for (auto particleIndex : *cell.getCellParticleIndices())
            {
                for (auto &neighbouringCellIndex : neighbours)
                {
                    for (auto interactingParticleIndex : *cells[neighbouringCellIndex].getCellParticleIndices())
                    {
                        // current particle is always within domain
                        Particle *currentParticle = &activeParticles[particleIndex];

                        // interacting particle can also be in halo
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
        if (!calcX)
            continue;

        // check if particle moved out of cell
        int cellIndex = computeCellIdx(particle);

        std::array<int, 3> idx3D = PContainer::convert1DTo3D(cellIndex, numCells);
        bool outOfBounds = idx3D[0] < 0 || idx3D[0] >= numCells[0] || idx3D[1] < 0 || idx3D[1] >= numCells[1] || idx3D[2] < 0 || idx3D[2] >= numCells[2];

        // particle out of bounds (not in domain or halo cell layer)
        if (outOfBounds)
        {
            particle.setHalo(true);
            cellRebuild = true;
            _simulationLogger->debug("Particle way out of bounds (" + std::to_string(cellIndex) + "): " + particle.toString());
        }
        else if (cellIndex != particle.getCellIdx())
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

const void LinkedCellParticleContainer::addParticle(std::array<double, 3> &x, std::array<double, 3> &v, double &m, double &epsilon, double &sigma, int &type)
{
    if (x[0] >= 0 && x[0] < domain[0] && x[1] >= 0 && x[1] < domain[1] && x[2] >= 0 && x[2] < domain[2])
    {

        activeParticles.emplace_back(x, v, m, epsilon, sigma, type);
        rebuildCells();
    }
}

const void LinkedCellParticleContainer::addParticle(std::array<double, 3> &x, std::array<double, 3> &v, double &m, double &epsilon, double &sigma, int &type, double &stiffness, double &averageBondLength)
{
    if (x[0] >= 0 && x[0] < domain[0] && x[1] >= 0 && x[1] < domain[1] && x[2] >= 0 && x[2] < domain[2])
    {

        activeParticles.emplace_back(x, v, m, epsilon, sigma, type, stiffness, averageBondLength);
        rebuildCells();
    }
}

const void LinkedCellParticleContainer::addParticle(std::array<double, 3> &x, std::array<double, 3> &v, double &m, double &epsilon, double &sigma, int &type, double &stiffness, double &averageBondLength, bool &fixed)
{
    if (x[0] >= 0 && x[0] < domain[0] && x[1] >= 0 && x[1] < domain[1] && x[2] >= 0 && x[2] < domain[2])
    {

        activeParticles.emplace_back(x, v, m, epsilon, sigma, type, stiffness, averageBondLength, fixed);
        rebuildCells();
    }
}


const void LinkedCellParticleContainer::addParticle(std::array<double, 3> &x, std::array<double, 3> &v, double &m, double &epsilon, double &sigma, bool &fixed, int &type)
{
    if (x[0] >= 0 && x[0] < domain[0] && x[1] >= 0 && x[1] < domain[1] && x[2] >= 0 && x[2] < domain[2])
    {

        activeParticles.emplace_back(x, v, m, epsilon, sigma, fixed, type);
        rebuildCells();
    }
}


const void LinkedCellParticleContainer::addParticle(std::array<double, 3> &x, std::array<double, 3> &v, std::array<double, 3> &f, std::array<double, 3> &old_f, double &m, double &epsilon, double &sigma, int &type, double &stiffness, double &averageBondLength)
{
    if (x[0] >= 0 && x[0] < domain[0] && x[1] >= 0 && x[1] < domain[1] && x[2] >= 0 && x[2] < domain[2])
    {

        activeParticles.emplace_back(x, v, f, old_f, m, epsilon, sigma, type, stiffness, averageBondLength);
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
            for (auto particleIndex : *cells[cellIndex].getCellParticleIndices())
            {
                Particle &p = activeParticles[particleIndex];
                double reflectingDistance = std::pow(2, 1.0 / 6) * p.getSigma() / 2;
                counterParticle.setEpsilon(p.getEpsilon());
                counterParticle.setSigma(p.getSigma());
                switch (b)
                {
                // left boundary, x-coordinate is 0
                case 0:
                    if (p.getX()[0] < reflectingDistance)
                    {
                        std::array<double, 3> counterX = {(-1) * p.getX()[0], p.getX()[1], p.getX()[2]};
                        counterParticle.setX(counterX);
                        f(p, counterParticle);
                    }
                    break;
                // right boundary, x-coordinate is domain size in x-direction
                case 1:
                    if (domain[0] - p.getX()[0] < reflectingDistance)
                    {
                        std::array<double, 3> counterX = {domain[0] + domain[0] - p.getX()[0], p.getX()[1], p.getX()[2]};
                        counterParticle.setX(counterX);
                        f(p, counterParticle);
                    }
                    break;
                // lower boundary, y-coordinate is 0
                case 2:
                    if (p.getX()[1] <= reflectingDistance)
                    {
                        std::array<double, 3> counterX = {p.getX()[0], (-1) * p.getX()[1], p.getX()[2]};
                        counterParticle.setX(counterX);
                        f(p, counterParticle);
                    }
                    break;
                // upper boundary, y-coordinate is domain size in y-direction
                case 3:
                    if (domain[1] - p.getX()[1] < reflectingDistance)
                    {
                        std::array<double, 3> counterX = {p.getX()[0], domain[1] + domain[1] - p.getX()[1], p.getX()[2]};
                        counterParticle.setX(counterX);
                        f(p, counterParticle);
                    }
                    break;
                // front boundary, z-coordinate is 0
                case 4:
                    if (p.getX()[2] <= reflectingDistance)
                    {
                        std::array<double, 3> counterX = {p.getX()[0], p.getX()[1], (-1) * p.getX()[2]};
                        counterParticle.setX(counterX);
                        f(p, counterParticle);
                    }
                    break;
                // back boundary, z-coordinate is domain size in z-direction
                case 5:
                    if (domain[2] - p.getX()[2] < reflectingDistance)
                    {
                        std::array<double, 3> counterX = {p.getX()[0], p.getX()[1], domain[2] + domain[2] - p.getX()[2]};
                        counterParticle.setX(counterX);
                        f(p, counterParticle);
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
                if (activeParticles[particleIndex].getFixed())
                    continue;
                Particle &toMirror = activeParticles[particleIndex];
                std::array<double, 3> newX = toMirror.getX() + mirroringOffset;
                int type = 0; 
                haloParticles.emplace_back(newX, toMirror.getV(), toMirror.getM(), toMirror.getEpsilon(), toMirror.getSigma(), type);
                Particle &mirrored = haloParticles.back();
                int cellIdx = computeCellIdx(mirrored);
                mirrored.setCellIdx(cellIdx);
                cells[cellIdx].insertParticleIndex(haloParticles.size() - 1);
                _simulationLogger->debug("Mirrored particle at " + std::to_string(newX[0]) + ", " + std::to_string(newX[1]) + ", " + std::to_string(newX[2]));
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