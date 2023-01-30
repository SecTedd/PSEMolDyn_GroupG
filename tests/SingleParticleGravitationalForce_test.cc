#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <vector>
#include "../src/model/DirectSumParticleContainer.h"
#include "../src/utils/ArrayUtils.h"
#include "../src/simulation/SingleParticleGravitationalForce.h"

// check correctness of Lennard-Jones-Forcecalculation against hand-calculated values
TEST(SingleParticleGravitationalForce, GravitationalForce)
{
    DirectSumParticleContainer pc = DirectSumParticleContainer();
    std::array<double, 3> x1 = {1, 0, 0};
    std::array<double, 3> v = {0, 0, 0};
    double m = 1;
    double sigma = 1;
    double epsilon = 5;
    int type = 1; 
    std::array<double, 3> g_grav = {0, -12, 0};
    pc.addParticle(x1, v, m, epsilon, sigma, type);

    // calculating new forces according to Lennard-Jones potential with hardcoded values epsilon=5 and sigma=1
    std::unique_ptr<SingleParticleGravitationalForce> calculation = std::make_unique<SingleParticleGravitationalForce>(g_grav);

    calculation->calculateForce(pc, 0);

    std::vector<Particle> &particles = pc.getActiveParticles();

    EXPECT_THAT(particles[0].getF(), testing::ElementsAre(0, -12, 0));
}