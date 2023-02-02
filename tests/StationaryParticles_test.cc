#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <memory>

#include "../src/model/Particle.h"
#include "../src/model/LinkedCellParticleContainer.h"
#include "../src/utils/ArrayUtils.h"

// Test that stationary particles dont experience any force
TEST(Particles, Stationary)
{
    BoundaryCondition o = BoundaryCondition::Outflow;
    BoundaryCondition r = BoundaryCondition::Reflecting;
    std::array<double, 3> domain = {80.0, 80.0, 1.0};
    std::array<BoundaryCondition, 6> bound = {r, r, r, r, o, o};
    double cutoff = 3;
    auto pc = std::make_shared<LinkedCellParticleContainer>(cutoff, domain, bound);

    std::array<double, 3> x1 = {0.0, 0.0, 0.0};
    std::array<double, 3> x2 = {0.8, 0.0, 0.0};
    std::array<double, 3> x3 = {0.7, 1.0, 0.0};
    std::array<double, 3> x4 = {0.0, 0.5, 0.0};
    std::array<double, 3> v = {0.0, 0.0, 0.0};
    double m = 1;
    double sigma = 1;
    double epsilon = 5;
    bool fixed = true;
    int type = 1;

    pc->addParticle(x1, v, m, epsilon, sigma, fixed, type);
    pc->addParticle(x2, v, m, epsilon, sigma, fixed, type);
    pc->addParticle(x3, v, m, epsilon, sigma, fixed, type);
    pc->addParticle(x4, v, m, epsilon, sigma, fixed, type);

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

    pc->iterateParticleInteractions(forceCalculationIteration);

    // validate that no force was apllied to particles
    for (auto &p: pc->getActiveParticles()) {
        EXPECT_TRUE(p.getF().at(0) == 0.0 && p.getF().at(1) == 0.0 && p.getF().at(2) == 0.0);
    }
}