#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "../src/model/LinkedCellParticleContainer.h"
#include "../src/utils/ArrayUtils.h"


/**
 * Test correct initialization of cells groups if compiled without OpenMP (no parallelization)
*/
TEST(CellGrouping, NoParallelization) {
    double cutoff = 3;
    std::array<double, 3> domain = {6, 6, 6};
    BoundaryCondition out = BoundaryCondition::Outflow;
    std::array<BoundaryCondition, 6> boundaries = {out, out, out, out, out, out};
    int parallelization = 0; //since tests are compiled without OpenMP this would automatically reset to 0 anyways
    LinkedCellParticleContainer pc = LinkedCellParticleContainer(cutoff, domain, boundaries, parallelization);
    std::vector<std::vector<int>> &cellGroups = pc.getCellGroups();
    EXPECT_EQ(pc.getParallel(), 0);
    EXPECT_EQ(cellGroups.size(), 1);
    EXPECT_EQ(cellGroups[0].size(), 8);
    EXPECT_THAT(cellGroups[0], testing::ElementsAre(21, 22, 25, 26, 37, 38, 41, 42));
}

/**
 * Test correct initalization of cells groups for parallelization strategy 1
 * Since tests are compiled without OpenMP, manually reset strategy to 1
*/
TEST(CellGrouping, Strategy1) {
    double cutoff = 3;
    std::array<double, 3> domain = {9, 9, 9};
    BoundaryCondition out = BoundaryCondition::Outflow;
    std::array<BoundaryCondition, 6> boundaries = {out, out, out, out, out, out};
    int parallelization = 0; //since tests are compiled without OpenMP this would automatically reset to 0 anyways
    LinkedCellParticleContainer pc = LinkedCellParticleContainer(cutoff, domain, boundaries, parallelization);
    //reset parallelization strategy to 1 and reinitialize cells for correct grouping
    pc.setParallel(1);
    pc.initializeCells(boundaries);
    std::vector<std::vector<int>> &cellGroups = pc.getCellGroups();
    EXPECT_EQ(pc.getParallel(), 1);
    EXPECT_EQ(cellGroups.size(), 18);
    for (auto &group : cellGroups) {
        EXPECT_EQ(group.capacity(), 2);
    }
    EXPECT_THAT(cellGroups[0], testing::ElementsAre(31, 81));
    EXPECT_THAT(cellGroups[1], testing::ElementsAre(32, 82));
    EXPECT_THAT(cellGroups[2], testing::ElementsAre(33, 83));
    EXPECT_THAT(cellGroups[3], testing::ElementsAre(36, 86));
    EXPECT_THAT(cellGroups[4], testing::ElementsAre(37, 87));
    EXPECT_THAT(cellGroups[5], testing::ElementsAre(38, 88));
    EXPECT_THAT(cellGroups[6], testing::ElementsAre(41, 91));
    EXPECT_THAT(cellGroups[7], testing::ElementsAre(42, 92));
    EXPECT_THAT(cellGroups[8], testing::ElementsAre(43, 93));
    EXPECT_THAT(cellGroups[9], testing::ElementsAre(56));
    EXPECT_THAT(cellGroups[10], testing::ElementsAre(57));
    EXPECT_THAT(cellGroups[11], testing::ElementsAre(58));
    EXPECT_THAT(cellGroups[12], testing::ElementsAre(61));
    EXPECT_THAT(cellGroups[13], testing::ElementsAre(62));
    EXPECT_THAT(cellGroups[14], testing::ElementsAre(63));
    EXPECT_THAT(cellGroups[15], testing::ElementsAre(66));
    EXPECT_THAT(cellGroups[16], testing::ElementsAre(67));
    EXPECT_THAT(cellGroups[17], testing::ElementsAre(68));
}

/**
 * Test correct initalization of cells groups for parallelization strategy 2
 * Since tests are compiled without OpenMP, manually reset strategy to 2
*/
TEST(CellGrouping, Strategy2) {
    double cutoff = 3;
    std::array<double, 3> domain = {12, 12, 12};
    BoundaryCondition out = BoundaryCondition::Outflow;
    std::array<BoundaryCondition, 6> boundaries = {out, out, out, out, out, out};
    int parallelization = 0; //since tests are compiled without OpenMP this would automatically reset to 0 anyways
    LinkedCellParticleContainer pc = LinkedCellParticleContainer(cutoff, domain, boundaries, parallelization);
    //reset parallelization strategy to 2 and reinitialize cells for correct grouping
    pc.setParallel(2);
    pc.initializeCells(boundaries);
    std::vector<std::vector<int>> &cellGroups = pc.getCellGroups();
    EXPECT_EQ(pc.getParallel(), 2);
    EXPECT_EQ(cellGroups.size(), 18);
    for (auto &group : cellGroups) {
        EXPECT_EQ(group.capacity(), 1);
    }
    //groups contain supercell indices
    EXPECT_THAT(cellGroups[0], testing::ElementsAre(0));
    EXPECT_THAT(cellGroups[1], testing::ElementsAre(1));
    EXPECT_THAT(cellGroups[2], testing::ElementsAre());
    EXPECT_THAT(cellGroups[3], testing::ElementsAre(2));
    EXPECT_THAT(cellGroups[4], testing::ElementsAre(3));
    EXPECT_THAT(cellGroups[5], testing::ElementsAre());
    EXPECT_THAT(cellGroups[6], testing::ElementsAre());
    EXPECT_THAT(cellGroups[7], testing::ElementsAre());
    EXPECT_THAT(cellGroups[8], testing::ElementsAre());
    EXPECT_THAT(cellGroups[9], testing::ElementsAre(4));
    EXPECT_THAT(cellGroups[10], testing::ElementsAre(5));
    EXPECT_THAT(cellGroups[11], testing::ElementsAre());
    EXPECT_THAT(cellGroups[12], testing::ElementsAre(6));
    EXPECT_THAT(cellGroups[13], testing::ElementsAre(7));
    EXPECT_THAT(cellGroups[14], testing::ElementsAre());
    EXPECT_THAT(cellGroups[15], testing::ElementsAre());
    EXPECT_THAT(cellGroups[16], testing::ElementsAre());
    EXPECT_THAT(cellGroups[17], testing::ElementsAre());

    std::vector<std::vector<int>> &superCells = pc.getSuperCells();
    EXPECT_EQ(superCells.size(), 8);
    for (auto &superCell : superCells) {
        EXPECT_EQ(superCell.capacity(), 8);
    }

    EXPECT_THAT(superCells[0], testing::ElementsAre(43, 44, 49, 50, 79, 80, 85, 86));
    EXPECT_THAT(superCells[1], testing::ElementsAre(45, 46, 51, 52, 81, 82, 87, 88));
    EXPECT_THAT(superCells[2], testing::ElementsAre(55, 56, 61, 62, 91, 92, 97, 98));
    EXPECT_THAT(superCells[3], testing::ElementsAre(57, 58, 63, 64, 93, 94, 99, 100));
    EXPECT_THAT(superCells[4], testing::ElementsAre(115, 116, 121, 122, 151, 152, 157, 158));
    EXPECT_THAT(superCells[5], testing::ElementsAre(117, 118, 123, 124, 153, 154, 159, 160));
    EXPECT_THAT(superCells[6], testing::ElementsAre(127, 128, 133, 134, 163, 164, 169, 170));
    EXPECT_THAT(superCells[7], testing::ElementsAre(129, 130, 135, 136, 165, 166, 171, 172));
}

/**
 * Test correct initalization of cells groups for parallelization strategy 3
 * Since tests are compiled without OpenMP, manually reset strategy to 3
*/
TEST(CellGrouping, Strategy3) {
    double cutoff = 3;
    std::array<double, 3> domain = {9, 9, 9};
    BoundaryCondition out = BoundaryCondition::Outflow;
    std::array<BoundaryCondition, 6> boundaries = {out, out, out, out, out, out};
    int parallelization = 0; //since tests are compiled without OpenMP this would automatically reset to 0 anyways
    LinkedCellParticleContainer pc = LinkedCellParticleContainer(cutoff, domain, boundaries, parallelization);
    //reset parallelization strategy to 3 and reinitialize cells for correct grouping
    pc.setParallel(3);
    pc.initializeCells(boundaries);

    std::vector<std::vector<int>> &cellGroups = pc.getCellGroups();
    EXPECT_EQ(pc.getParallel(), 3);
    EXPECT_EQ(cellGroups.size(), 1);
    EXPECT_EQ(cellGroups[0].size(), 8);

    std::vector<std::vector<int>> &superCells = pc.getSuperCells();
    EXPECT_EQ(superCells.size(), 8);
    for (auto &superCell : superCells) {
        EXPECT_EQ(superCell.capacity(), 8);
    }
    
    EXPECT_THAT(superCells[0], testing::ElementsAre(31, 32, 36, 37, 56, 57, 61, 62));
    EXPECT_THAT(superCells[1], testing::ElementsAre(33, 38, 58, 63));
    EXPECT_THAT(superCells[2], testing::ElementsAre(41, 42, 66, 67));
    EXPECT_THAT(superCells[3], testing::ElementsAre(43, 68));
    EXPECT_THAT(superCells[4], testing::ElementsAre(81, 82, 86, 87));
    EXPECT_THAT(superCells[5], testing::ElementsAre(83, 88));
    EXPECT_THAT(superCells[6], testing::ElementsAre(91, 92));
    EXPECT_THAT(superCells[7], testing::ElementsAre(93));
    
}

/**
 * Test correct force interaction between particles for groups corresponding to parallelization strategy 1
 * NOT test of actual parallel execution
 * Since tests are compiled without OpenMP, manually reset strategy to 1
*/
TEST(ParallelParticleInteractions, Strategy1) {
    double cutoff = 3;
    std::array<double, 3> domain = {9, 9, 9};
    BoundaryCondition out = BoundaryCondition::Outflow;
    std::array<BoundaryCondition, 6> boundaries = {BoundaryCondition::Reflecting, out, out, out, out, out};
    int parallelization = 0; //since tests are compiled without OpenMP this would automatically reset to 0 anyways
    LinkedCellParticleContainer pc = LinkedCellParticleContainer(cutoff, domain, boundaries, parallelization);
    //reset parallelization strategy to 1 and reinitialize cells for correct grouping
    pc.setParallel(1);
    pc.initializeCells(boundaries);

    std::array<double, 3> x1 = {0.5, 0, 0}; // only interaction with border
    std::array<double, 3> x2 = {4, 2.5, 0}; // only interactinon with x3
    std::array<double, 3> x3 = {4, 3.5, 0}; // only interaction with x2
    std::array<double, 3> x4 = {8, 8, 2.5}; // only interaction with x5
    std::array<double, 3> x5 = {8, 8, 3.5}; // only interaction with x6

    std::array<double, 3> v = {0, 0, 0};
    double m = 1;
    double epsilon = 1;
    double sigma = 1;

    pc.reserveMemoryForParticles(5);
    pc.addParticle(x1, v, m, epsilon, sigma);
    pc.addParticle(x2, v, m, epsilon, sigma);
    pc.addParticle(x3, v, m, epsilon, sigma);
    pc.addParticle(x4, v, m, epsilon, sigma);
    pc.addParticle(x5, v, m, epsilon, sigma);

    // force calculation function
    std::function<void(Particle &, Particle &)> forceCalculationIteration = [](Particle &p1, Particle &p2)
    {
        double distance = ArrayUtils::L2Norm(p1.getX() - p2.getX());

        // Reduce number of operation by reusing previous results
        double pow1 = 1 / distance;
        double pow2 = pow1 * pow1;
        double pow4 = pow2 * pow2;
        double pow6 = pow4 * pow2;
        double pow12 = pow6 * pow6;

        // Lennard-Jones force
        std::array<double, 3> f_ij = (-24 * 5 / pow(distance, 2)) * (pow6 - 2 * pow12) * (p1.getX() - p2.getX());
        std::array<double, 3> f_ji = -1 * f_ij;

        p1.addF(f_ij);
        p2.addF(f_ji);
    };

    pc.iterateParticleInteractions(forceCalculationIteration);

    std::vector<Particle> &particles = pc.getActiveParticles();

    EXPECT_EQ(pc.getParallel(), 1);

    EXPECT_THAT(particles[0].getF(), testing::ElementsAre(120, 0, 0));
    EXPECT_THAT(particles[1].getF(), testing::ElementsAre(0, -120, 0));
    EXPECT_THAT(particles[2].getF(), testing::ElementsAre(0, 120, 0));
    EXPECT_THAT(particles[3].getF(), testing::ElementsAre(0, 0, -120));
    EXPECT_THAT(particles[4].getF(), testing::ElementsAre(0, 0, 120));
}

/**
 * Test correct force interaction between particles for groups corresponding to parallelization strategy 2
 * NOT test of actual parallel execution
 * Since tests are compiled without OpenMP, manually reset strategy to 1
*/
TEST(ParallelParticleInteractions, Strategy2) {
    double cutoff = 3;
    std::array<double, 3> domain = {9, 9, 9};
    BoundaryCondition out = BoundaryCondition::Outflow;
    std::array<BoundaryCondition, 6> boundaries = {BoundaryCondition::Reflecting, out, out, out, out, out};
    int parallelization = 0; //since tests are compiled without OpenMP this would automatically reset to 0 anyways
    LinkedCellParticleContainer pc = LinkedCellParticleContainer(cutoff, domain, boundaries, parallelization);
    //reset parallelization strategy to 1 and reinitialize cells for correct grouping
    pc.setParallel(2);
    pc.initializeCells(boundaries);

    std::array<double, 3> x1 = {0.5, 0, 0}; // only interaction with border
    std::array<double, 3> x2 = {4, 2.5, 0}; // only interactinon with x3
    std::array<double, 3> x3 = {4, 3.5, 0}; // only interaction with x2
    std::array<double, 3> x4 = {8, 8, 2.5}; // only interaction with x5
    std::array<double, 3> x5 = {8, 8, 3.5}; // only interaction with x6

    std::array<double, 3> v = {0, 0, 0};
    double m = 1;
    double epsilon = 1;
    double sigma = 1;

    pc.reserveMemoryForParticles(5);
    pc.addParticle(x1, v, m, epsilon, sigma);
    pc.addParticle(x2, v, m, epsilon, sigma);
    pc.addParticle(x3, v, m, epsilon, sigma);
    pc.addParticle(x4, v, m, epsilon, sigma);
    pc.addParticle(x5, v, m, epsilon, sigma);

    // force calculation function
    std::function<void(Particle &, Particle &)> forceCalculationIteration = [](Particle &p1, Particle &p2)
    {
        double distance = ArrayUtils::L2Norm(p1.getX() - p2.getX());

        // Reduce number of operation by reusing previous results
        double pow1 = 1 / distance;
        double pow2 = pow1 * pow1;
        double pow4 = pow2 * pow2;
        double pow6 = pow4 * pow2;
        double pow12 = pow6 * pow6;

        // Lennard-Jones force
        std::array<double, 3> f_ij = (-24 * 5 / pow(distance, 2)) * (pow6 - 2 * pow12) * (p1.getX() - p2.getX());
        std::array<double, 3> f_ji = -1 * f_ij;

        p1.addF(f_ij);
        p2.addF(f_ji);
    };

    pc.iterateParticleInteractions(forceCalculationIteration);

    std::vector<Particle> &particles = pc.getActiveParticles();

    EXPECT_EQ(pc.getParallel(), 2);

    EXPECT_THAT(particles[0].getF(), testing::ElementsAre(120, 0, 0));
    EXPECT_THAT(particles[1].getF(), testing::ElementsAre(0, -120, 0));
    EXPECT_THAT(particles[2].getF(), testing::ElementsAre(0, 120, 0));
    EXPECT_THAT(particles[3].getF(), testing::ElementsAre(0, 0, -120));
    EXPECT_THAT(particles[4].getF(), testing::ElementsAre(0, 0, 120));
}

/**
 * No unit test for interactions of parallel strategy 3
 * because tests are not compiled with OpenMP, omp tasks cannot be executed
*/
