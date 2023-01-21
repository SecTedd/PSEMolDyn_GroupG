#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <vector>
#include "../src/model/Cuboid.h"
#include "../src/model/ProgramParameters.h"
#include "../src/model/DirectSumParticleContainer.h"
#include "../src/utils/ArrayUtils.h"
#include "../src/utils/ParticleGenerator.h"
#include "../src/simulation/LennardJonesForceHarmonic.h"
#include "../src/simulation/TemporalSingleParticleForce.h"

// check correctness of Lennard-Jones-Force with the harmonic potential calculation against hand-calculated values
TEST(Membrane, LennardJonesForceHarmonicRepulsive)
{
    DirectSumParticleContainer pc = DirectSumParticleContainer();
    std::array<double, 3> x1 = {1, 1, 0};
    std::array<double, 3> x2 = {2, 1, 0};
    std::array<double, 3> v = {0, 0, 0};
    // distance of particles smaller than average bond length, particles repell each other for both harmonic potential and lennard jones
    double h = 1;
    double m = 1;
    double sigma = 1;
    double epsilon = 5;
    int type = 0;
    double stiffness = 300;
    double averageBondLength = 2;
    auto p1 = Particle(x1, v, m, epsilon, sigma, type, stiffness, averageBondLength);
    auto p2 = Particle(x2, v, m, epsilon, sigma, type, stiffness, averageBondLength);
    std::vector<int> parallel = std::vector<int>{1};
    p1.setParallelNeighbours(parallel);
    pc.getActiveParticles().emplace_back(p1);
    pc.getActiveParticles().emplace_back(p2);

    auto force = LennardJonesForceHarmonic();

    force.calculateForce(pc);

    EXPECT_THAT(pc.getActiveParticles()[0].getF(), testing::ElementsAre(-420, 0, 0));
    EXPECT_THAT(pc.getActiveParticles()[1].getF(), testing::ElementsAre(420, 0, 0));

    pc.resetParticles();

    std::array<double, 3> x4 = {2, 2, 0};
    auto p3 = Particle(x1, v, m, epsilon, sigma, type, stiffness, averageBondLength);
    auto p4 = Particle(x4, v, m, epsilon, sigma, type, stiffness, averageBondLength);
    std::vector<int> diagonal = std::vector<int>{1};
    p3.setDiagonalNeighbours(diagonal);
    pc.getActiveParticles().emplace_back(p3);
    pc.getActiveParticles().emplace_back(p4);

    force.calculateForce(pc);

    EXPECT_THAT(pc.getActiveParticles()[0].getF(), testing::ElementsAre(-300, -300, 0));
    EXPECT_THAT(pc.getActiveParticles()[1].getF(), testing::ElementsAre(300, 300, 0));
}

// check correctness of Lennard-Jones-Force with the harmonic potential calculation against hand-calculated values
TEST(Membrane, LennardJonesForceHarmonicAttracting)
{
    DirectSumParticleContainer pc = DirectSumParticleContainer();
    std::array<double, 3> x1 = {1, 1, 0};
    std::array<double, 3> x2 = {4, 1, 0};
    std::array<double, 3> v = {0, 0, 0};
    // distance of particles smaller than average bond length, particles repell each other for both harmonic potential and lennard jones
    double h = 1;
    double m = 1;
    double sigma = 1;
    double epsilon = 5;
    int type = 0;
    double stiffness = 300;
    double averageBondLength = 2;
    auto p1 = Particle(x1, v, m, epsilon, sigma, type, stiffness, averageBondLength);
    auto p2 = Particle(x2, v, m, epsilon, sigma, type, stiffness, averageBondLength);
    std::vector<int> parallel = std::vector<int>{1};
    p1.setParallelNeighbours(parallel);
    pc.getActiveParticles().emplace_back(p1);
    pc.getActiveParticles().emplace_back(p2);

    auto force = LennardJonesForceHarmonic();

    force.calculateForce(pc);

    EXPECT_THAT(pc.getActiveParticles()[0].getF(), testing::ElementsAre(300, 0, 0));
    EXPECT_THAT(pc.getActiveParticles()[1].getF(), testing::ElementsAre(-300, 0, 0));

    pc.resetParticles();

    averageBondLength = 1;
    std::array<double, 3> x4 = {3, 3, 0};
    auto p3 = Particle(x1, v, m, epsilon, sigma, type, stiffness, averageBondLength);
    auto p4 = Particle(x4, v, m, epsilon, sigma, type, stiffness, averageBondLength);
    std::vector<int> diagonal = std::vector<int>{1};
    p3.setDiagonalNeighbours(diagonal);
    pc.getActiveParticles().emplace_back(p3);
    pc.getActiveParticles().emplace_back(p4);

    force.calculateForce(pc);

    EXPECT_THAT(pc.getActiveParticles()[0].getF(), testing::ElementsAre(300, 300, 0));
    EXPECT_THAT(pc.getActiveParticles()[1].getF(), testing::ElementsAre(-300, -300, 0));
}

TEST(Membrane, ParticleNeighbours)
{
    DirectSumParticleContainer pc = DirectSumParticleContainer();
    std::array<double, 3> x = {1, 1, 0};
    std::array<double, 3> v = {0, 0, 0};
    std::array<int, 3> n = {3, 3, 1};
    double h = 2.2;
    double m = 1;
    double sigma = 1;
    double epsilon = 5;
    int type = 0;
    auto cuboid = Cuboid(x, n, h, m, v, epsilon, sigma, type);

    ParticleGenerator::generateCuboid(pc, cuboid, true);

    EXPECT_THAT(pc.getActiveParticles()[0].getParallelNeighbours(), testing::ElementsAre(3, 1));
    EXPECT_THAT(pc.getActiveParticles()[0].getDiagonalNeighbours(), testing::ElementsAre(4));
    EXPECT_THAT(pc.getActiveParticles()[4].getParallelNeighbours(), testing::ElementsAre(7, 5));
    EXPECT_THAT(pc.getActiveParticles()[4].getDiagonalNeighbours(), testing::ElementsAre(6, 8));
    EXPECT_THAT(pc.getActiveParticles()[8].getParallelNeighbours(), testing::ElementsAre());
    EXPECT_THAT(pc.getActiveParticles()[8].getDiagonalNeighbours(), testing::ElementsAre());
}

TEST(Membrane, TemporalSingleParticleForce)
{
    DirectSumParticleContainer pc = DirectSumParticleContainer();
    std::array<double, 3> x1 = {1, 1, 1};
    std::array<double, 3> x2 = {2, 2, 1};
    std::array<double, 3> v = {0, 0, 0};
    double m = 1;
    double sigma = 1;
    double epsilon = 5;
    int type = 1;
    std::array<double, 3> f = {1, 1, 1};
    std::vector<int> particles = {0};
    double endTime = 1; 
    pc.addParticle(x1, v, m, epsilon, sigma, type);
    pc.addParticle(x2, v, m, epsilon, sigma, type);

    auto force = TemporalSingleParticleForce(f, endTime, particles); 

    force.calculateForce(pc, 2); 
    EXPECT_THAT(pc.getActiveParticles()[0].getF(), testing::ElementsAre(0, 0, 0));
    EXPECT_THAT(pc.getActiveParticles()[1].getF(), testing::ElementsAre(0, 0, 0));


    force.calculateForce(pc, 0.5); 
    EXPECT_THAT(pc.getActiveParticles()[0].getF(), testing::ElementsAre(1, 1, 1));
    EXPECT_THAT(pc.getActiveParticles()[1].getF(), testing::ElementsAre(0, 0, 0));
}