#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "../src/simulation/DVProfileCalculator.h"
#include "../src/model/LinkedCellParticleContainer.h"

// basic test for particle bin calculation
TEST(DVCalculator, ParticlesPerBin) 
{
    BoundaryCondition o = BoundaryCondition::Outflow;
    std::array<double, 3> domain = {10.0, 10.0, 1.0};
    std::array<BoundaryCondition, 6> bound = {o, o, o, o, o, o};
    double cutoff = 3;
    auto pc = std::make_shared<LinkedCellParticleContainer>(cutoff, domain, bound, 0);

    std::array<double, 3> x1 = {2.0, 0.0, 0.0};
    std::array<double, 3> x2 = {8.0, 0.0, 0.0};
    std::array<double, 3> v = {100.0, 100.0, 0.0};
    double m = 1;
    double sigma = 1;
    double epsilon = 5;
    int type = 1;
    
    pc->addParticle(x1, v, m, epsilon, sigma, type);
    pc->addParticle(x2, v, m, epsilon, sigma, type);

    auto dv_calc = std::make_shared<DVProfileCalculator>(pc, 2, domain);
    auto result = dv_calc->calculate();

    EXPECT_TRUE(result.size() == 2);
    EXPECT_TRUE(result.at(0) == 1 && result.at(1) == 1);
}

// tests the correct value of the internal avg attribute
TEST(DVCalculator, AverageParticlesPerBin)
{
    BoundaryCondition o = BoundaryCondition::Outflow;
    std::array<double, 3> domain = {10.0, 10.0, 1.0};
    std::array<BoundaryCondition, 6> bound = {o, o, o, o, o, o};
    double cutoff = 3;
    auto pc = std::make_shared<LinkedCellParticleContainer>(cutoff, domain, bound, 0);

    std::array<double, 3> x1 = {2.0, 0.0, 0.0};
    std::array<double, 3> x2 = {8.0, 0.0, 0.0};
    std::array<double, 3> v = {100.0, 100.0, 0.0};
    double m = 1;
    double sigma = 1;
    double epsilon = 5;
    int type = 1;

    pc->addParticle(x1, v, m, epsilon, sigma, type);
    pc->addParticle(x2, v, m, epsilon, sigma, type);

    auto dv_calc = std::make_shared<DVProfileCalculator>(pc, 2, domain);
    auto result = dv_calc->calculate();
    result = dv_calc->calculate();
    result = dv_calc->calculate();
    result = dv_calc->calculate();
    result = dv_calc->calculate();
    result = dv_calc->calculate();

    EXPECT_TRUE(dv_calc->getAvg().size() == 2);
    EXPECT_TRUE(dv_calc->getAvg().at(0) == 6 && dv_calc->getAvg().at(1) == 6);
}