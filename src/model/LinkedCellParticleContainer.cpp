
#include "LinkedCellParticleContainer.h"
#include "../utils/PContainer.h"
#include "../utils/ArrayUtils.h"

#include <iostream>
#include <cmath>
#include <omp.h>

LinkedCellParticleContainer::LinkedCellParticleContainer(double cutoff, std::array<double, 3> &domain, std::array<BoundaryCondition, 6> &domainBoundaries, int parallel)
{
    _memoryLogger = spdlog::get("memory_logger");
    _memoryLogger->info("LinkedCellParticleContainer generated!");
    _simulationLogger = spdlog::get("simulation_logger");

    this->domain = domain;
    this->cutoff = cutoff;

// use given parallelization strategy if compiled with OpenMP
// no parallelization otherwise
    #ifdef _OPENMP
    this->parallel = parallel;
    #else
    this->parallel = 0;
    #endif

    _memoryLogger->warn("Parallel set to " + std::to_string(this->parallel));

    initializeCells(domainBoundaries);
}

LinkedCellParticleContainer::~LinkedCellParticleContainer()
{
    _memoryLogger->info("LinkedCellParticleContainer destructed!");
}

const void LinkedCellParticleContainer::initializeCells(std::array<BoundaryCondition, 6> &domainBoundaries)
{
    //make sure to start from empty vectors
    cells.clear();
    cellGroups.clear();
    superCells.clear();

    // smallest cellsize bigger than cutoff plus boundaries
    int numberOfXCells = static_cast<int>(std::floor(domain[0] / cutoff));
    int numberOfYCells = static_cast<int>(std::floor(domain[1] / cutoff));
    int numberOfZCells = static_cast<int>(std::floor(domain[2] / cutoff));

    // for 2D we need exactly one layer, which is not guaranteed above
    if (numberOfZCells == 0)
        numberOfZCells = 1;

    // add surrounding halo cells
    numCells = {numberOfXCells + 2, numberOfYCells + 2, numberOfZCells + 2};

    //if parallel = 2 use clustered supercells
    if (parallel == 2 || parallel == 3)
        numInteractingCells = {(numberOfXCells + 1) / 2, (numberOfYCells + 1) / 2, (numberOfZCells + 1) / 2};
    else
        numInteractingCells = {numberOfXCells, numberOfYCells, numberOfZCells};

    double sizeX = domain[0] / numberOfXCells;
    double sizeY = domain[1] / numberOfYCells;
    double sizeZ = domain[2] / numberOfZCells;

    // cell size has to be at least 1, otherwise we get problems with division by 0
    if (sizeZ == 0)
        sizeZ = 1;

    cellSize = {sizeX, sizeY, sizeZ};

    _memoryLogger->warn("Cell size: " + std::to_string(sizeX) + ", " + std::to_string(sizeY) + ", " + std::to_string(sizeZ));
    _memoryLogger->warn("Num cells: " + std::to_string(numCells[0]) + ", " + std::to_string(numCells[1]) + ", " + std::to_string(numCells[2]));
    _memoryLogger->warn("Num supercells: " + std::to_string(numInteractingCells[0]) + ", " + std::to_string(numInteractingCells[1]) + ", " + std::to_string(numInteractingCells[2]));

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
    //Initialization needed for parallelism
    reserveGroups();
    initializeParallelGroups();
}

void LinkedCellParticleContainer::initializeParallelGroups() {
    for (unsigned int i = 0; i < cells.size(); i++) {
        if (cells[i].getType() == CellType::InnerCell || cells[i].getType() == CellType::BoundaryCell) {
            if (parallel == 2 || parallel == 3) {
                superCells[computeCellGroup(i)].emplace_back(i);
            }
            else {
                cellGroups[computeCellGroup(i)].emplace_back(i);
            }
        }
    }
    if (parallel == 2 || parallel == 3) {
        for (unsigned int i = 0; i < superCells.size(); i++) {
            cellGroups[computeSupercellGroup(i)].emplace_back(i);
        }
    }
}


void LinkedCellParticleContainer::reserveGroups()
{
    // default values if no parallelization is applied
    // code is not compiled with OpenMP so cells can be all in one group
    // if tasks are used parallelization doesn't need groups
    int numGroups = 1;
    int maxCellsPerGroup = numInteractingCells[0] * numInteractingCells[1] * numInteractingCells[2];

    // constant number of groups holding either supercells (parallel = 2) or cells directly (parallel = 1)
    if (parallel == 1 || parallel == 2)
    {
        // distinguish 2D and 3D cases
        numGroups = numInteractingCells[2] > 1 ? 18 : 6;
        int offsetX = numInteractingCells[0] % 3 > 0 ? 1 : 0;
        int offsetY = numInteractingCells[1] % 3 > 0 ? 1 : 0;
        maxCellsPerGroup = numInteractingCells[2] > 1
                               ? (numInteractingCells[0] / 3 + offsetX) * (numInteractingCells[1] / 3 + offsetY) * (numInteractingCells[2] / 2 + numInteractingCells[2] % 2)
                               : (numInteractingCells[0] / 3 + offsetX) * (numInteractingCells[1] / 2 + numInteractingCells[1] % 2);
    }

    // reserve memory for groups & cells in groups
    // initialize group vectors
    cellGroups.reserve(numGroups);
    for (int i = 0; i < numGroups; i++)
    {
        std::vector<int> group;
        cellGroups.push_back(group);
        cellGroups[i].reserve(maxCellsPerGroup);
    }

    if (parallel == 2 || parallel == 3) {
        int maxInnerCells = numInteractingCells[2] > 1 ? 8 : 4;
        int totalSuperCells = numInteractingCells[0] * numInteractingCells[1] * numInteractingCells[2];
        superCells.reserve(totalSuperCells);
        for (int i = 0; i < totalSuperCells; i++) {
            std::vector<int> superCell;
            superCells.push_back(superCell);
            superCells[i].reserve(maxInnerCells);
        }
    }
}

const int LinkedCellParticleContainer::computeSupercellGroup(int cellIdx) {
    //all supercells in one group
    if (parallel == 3) 
        return 0;


    int numGroupsX = 3;
    //distinguish 2D and 3D cas
    int numGroupsY = numInteractingCells[2] > 1 ? 3 : 2;
    int numGroupsZ = numInteractingCells[2] > 1 ? 2 : 1;
    std::array<int, 3> index3D = PContainer::convert1DTo3D(cellIdx, numInteractingCells);
    int groupX = index3D[0] % numGroupsX;
    int groupY = index3D[1] % numGroupsY;
    int groupZ = index3D[2] % numGroupsZ;

    return groupX + numGroupsX * groupY + numGroupsX * numGroupsY * groupZ;
}

const int LinkedCellParticleContainer::computeCellGroup(int cellIdx)
{
    int groupX = 0, groupY = 0, groupZ = 0;
    int numGroupsX = 0, numGroupsY = 0;

    if (parallel == 1)
    {
        numGroupsX = 3;

        // distinguish 2D and 3D case
        numGroupsY = (numCells[2] - 2) > 1 ? 3 : 2;
        int numGroupsZ = (numCells[2] - 2) > 1 ? 2 : 1;

        std::array<int, 3> index3D = PContainer::convert1DTo3D(cellIdx, numCells);
        groupX = (index3D[0] - 1) % numGroupsX;
        groupY = (index3D[1] - 1) % numGroupsY;
        groupZ = (index3D[2] - 1) % numGroupsZ;
    }

    else if (parallel == 2 || parallel == 3) {
        numGroupsX = (numCells[0] - 1) / 2;
        numGroupsY = (numCells[1] - 1) / 2;
        std::array<int, 3> index3D = PContainer::convert1DTo3D(cellIdx, numCells);

        groupX = (index3D[0] - 1) / 2;
        groupY = (index3D[1] - 1) / 2;
        groupZ = (index3D[2] - 1) / 2;
    }

    return groupX + numGroupsX * groupY + numGroupsX * numGroupsY * groupZ;
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

inline void LinkedCellParticleContainer::intraCellInteraction(int cellIndex, std::function<void(Particle &, Particle &)> f)
{
    std::vector<int> *particleIndices = cells[cellIndex].getCellParticleIndices();
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
}

inline void LinkedCellParticleContainer::interCellInteraction(int i, int j, std::function<void(Particle &, Particle &)> f)
{
    //std::cout << "Neighbour interaction between " << i << " and " << j << std::endl;
    // particle interactions between different cells
    for (auto particleIndex : *cells[i].getCellParticleIndices())
    {
        for (auto interactingParticleIndex : *cells[j].getCellParticleIndices())
        {
            // current particle is always within domain
            Particle *currentParticle = &activeParticles[particleIndex];

            // interacting particle can also be in halo
            Particle *interactingParticle;
            if (cells[j].getType() == CellType::HaloCell || cells[j].getType() == CellType::PeriodicHaloCell)
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

void LinkedCellParticleContainer::directCellInteraction(std::function<void(Particle &, Particle &)> f) {
    // interaction between different cells
    // if compiled without OpenMP no parallelization of cell interactions
    for (auto &group : cellGroups)
    {
        // interactions of cells in one group can be executed in parallel
        // no parallelization if not compiled with openMP
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif 
        for (int idx : group)
        {
            for (auto neighbour : cells[idx].getDomainNeighbours())
            {
                interCellInteraction(idx, neighbour, f);
            }
        }
    }
}

void LinkedCellParticleContainer::nestedCellInteraction(std::function<void(Particle &, Particle &)> f) {
    for (auto &group : cellGroups) {
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (int superCellIdx : group) {
            for (int cellIdx : superCells[superCellIdx]) {
                for (int neighbour : cells[cellIdx].getDomainNeighbours()) {
                    interCellInteraction(cellIdx, neighbour, f);
                }
            }
        }
    }
}

void LinkedCellParticleContainer::forkJoin(std::function<void(Particle &, Particle &)> f)
{
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (long unsigned int i = 0; i < cells.size(); i++)
    {
        if (cells[i].getType() == CellType::BoundaryCell || cells[i].getType() == CellType::InnerCell) {
            //interaction within cell
            intraCellInteraction(i, f);

            //interaction with halo neighbours 
            //race conditions in halo particle forces don't matter
            for (auto j : cells[i].getPeriodicHaloNeighbours())
                interCellInteraction(i, j, f);
        }
    }

    if (parallel == 2)
        nestedCellInteraction(f);
    else
        directCellInteraction(f);
}

void LinkedCellParticleContainer::taskModel(std::function<void(Particle &, Particle &)> f) {
    #ifdef _OPENMP
    #pragma omp parallel
    #pragma omp single
    {

    for (unsigned int superCellIdx = 0; superCellIdx < superCells.size(); superCellIdx++) {
        //one task per supercell
        std::vector<int> neighbouringSupercells = PContainer::getNeigbhourGroupsNewton(superCellIdx, numInteractingCells);
        
        int numDependencies = neighbouringSupercells.size();
        std::vector<int> *dependencies[numDependencies];
        for (int j = 0; j < numDependencies; j++) {
            dependencies[j] = &superCells[neighbouringSupercells[j]];
        }
        std::vector<int> *superCell = &superCells[superCellIdx];

        #pragma omp task depend(inout: superCell, dependencies[0:numDependencies]) 
        {
            for (auto cellIdx : *superCell) {
                ParticleCell *cell = &cells[cellIdx];

                intraCellInteraction(cellIdx, f);

                std::vector<int> neighbours = cell->getDomainNeighbours();
                if (cell->getType() == CellType::BoundaryCell)
                {
                    //boundary cells might have to interact with their halo neighbours
                    neighbours.insert(neighbours.end(), cell->getPeriodicHaloNeighbours().begin(), cell->getPeriodicHaloNeighbours().end());
                }
                for (auto j : neighbours) {
                    interCellInteraction(cellIdx, j, f);
                }
            }
        }
    }
    }
    #pragma omp taskwait
    #endif
}

const void LinkedCellParticleContainer::iterateParticleInteractions(std::function<void(Particle &, Particle &)> f)
{
    // first we need to clear the halo
    clearHalo();

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

    if (parallel == 3)
        taskModel(f);
    else
        forkJoin(f);
}

const void LinkedCellParticleContainer::iterateParticles(std::function<void(Particle &)> f, bool calcX)
{
    // no deletions in the active particle vector, meaning it suffices to delete invalid particles from the cells and add the particles elsewhere
    bool cellUpdate = false;
    // deletions in the active particle vector, meaning indices within the cells are rendered invalid
    bool cellRebuild = false;

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif

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

        activeParticles.emplace_back(x, v, m, epsilon, sigma, type);
        rebuildCells();
    }
}

const void LinkedCellParticleContainer::addParticle(std::array<double, 3> &x, std::array<double, 3> &v, std::array<double, 3> &f, std::array<double, 3> &old_f, double &m, double &epsilon, double &sigma, int &type)
{
    if (x[0] >= 0 && x[0] < domain[0] && x[1] >= 0 && x[1] < domain[1] && x[2] >= 0 && x[2] < domain[2])
    {

        activeParticles.emplace_back(x, v, f, old_f, m, epsilon, sigma, type);
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
                Particle &toMirror = activeParticles[particleIndex];
                std::array<double, 3> newX = toMirror.getX() + mirroringOffset;
                haloParticles.emplace_back(newX, toMirror.getV(), toMirror.getM(), toMirror.getEpsilon(), toMirror.getSigma());
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

std::vector<std::vector<int>> &LinkedCellParticleContainer::getCellGroups() { return cellGroups; }

std::vector<std::vector<int>> &LinkedCellParticleContainer::getSuperCells() { return superCells; }

void LinkedCellParticleContainer::setParallel(int parallel) { this->parallel = parallel; }

const int LinkedCellParticleContainer::getParallel() { return parallel; }

