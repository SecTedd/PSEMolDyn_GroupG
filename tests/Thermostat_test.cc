#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <memory>

#include "../src/simulation/Thermostat.h"
#include "../src/model/Particle.h"
#include "../src/model/LinkedCellParticleContainer.h"


// Basic test where temp is instantly set to 40
TEST(Thermostat, NoDeltaTo40) 
{
    BoundaryCondition o = BoundaryCondition::Outflow;
    BoundaryCondition r = BoundaryCondition::Reflecting;
    std::array<double, 3> domain = {80.0, 80.0, 1.0};
    std::array<BoundaryCondition, 6> bound = {r, r, r, r, o, o};
    double cutoff = 3;
    int parallelization = 0; //since tests are compiled without OpenMP this would automatically reset to 0 anyways
    auto pc = std::make_shared<LinkedCellParticleContainer>(cutoff, domain, bound, parallelization);

    std::array<double, 3> x = {0.0, 0.0, 0.0};
    std::array<double, 3> v = {100.0, 100.0, 100.0};
    double m = 1;
    double sigma = 1;
    double epsilon = 5;
    
    pc->addParticle(x, v, m, epsilon, sigma);

    Thermostat t = Thermostat(pc, 100.0, 2);
    t.setTargetTemperature(40.0);

    t.apply();
    
    EXPECT_TRUE(std::abs(t.calculateCurrentTemperature() - 40) < .5);
}

// Heating
TEST(Thermostat, Heating) 
{
    BoundaryCondition o = BoundaryCondition::Outflow;
    BoundaryCondition r = BoundaryCondition::Reflecting;
    std::array<double, 3> domain = {80.0, 80.0, 1.0};
    std::array<BoundaryCondition, 6> bound = {r, r, r, r, o, o};
    double cutoff = 3;
    int parallelization = 0; //since tests are compiled without OpenMP this would automatically reset to 0 anyways
    auto pc = std::make_shared<LinkedCellParticleContainer>(cutoff, domain, bound, parallelization);

    std::array<double, 3> x1 = {1.0, 1.0, 0.0};
    std::array<double, 3> x2 = {79.0, 79.0, 79.0};
    std::array<double, 3> v = {1.0, 1.0, 1.0};
    double m = 1;
    double sigma = 1;
    double epsilon = 5;
    
    pc->addParticle(x1, v, m, epsilon, sigma);
    pc->addParticle(x2, v, m, epsilon, sigma);

    Thermostat t = Thermostat(pc, 10.0, 2);
    t.setTargetTemperature(1000.0);
    t.setTemperatureDelta(0.5);

    for (int i = 0; i < 1000; i++) {
        double oldTemp = t.calculateCurrentTemperature();
        t.apply();
        double newTemp = t.calculateCurrentTemperature();
        EXPECT_TRUE(oldTemp <= newTemp) << "Old temp: " << oldTemp << " New temp: " << newTemp;
    }
}

// Cooling
TEST(Thermostat, Cooling) 
{
    BoundaryCondition o = BoundaryCondition::Outflow;
    BoundaryCondition r = BoundaryCondition::Reflecting;
    std::array<double, 3> domain = {80.0, 80.0, 1.0};
    std::array<BoundaryCondition, 6> bound = {r, r, r, r, o, o};
    double cutoff = 3;
    int parallelization = 0; //since tests are compiled without OpenMP this would automatically reset to 0 anyways
    auto pc = std::make_shared<LinkedCellParticleContainer>(cutoff, domain, bound, parallelization);

    std::array<double, 3> x1 = {1.0, 1.0, 0.0};
    std::array<double, 3> x2 = {79.0, 79.0, 79.0};
    std::array<double, 3> v = {1000.0, 1000.0, 1000.0};
    double m = 1;
    double sigma = 1;
    double epsilon = 5;
    
    pc->addParticle(x1, v, m, epsilon, sigma);
    pc->addParticle(x2, v, m, epsilon, sigma);

    Thermostat t = Thermostat(pc, 1000.0, 2);
    t.setTargetTemperature(10.0);
    t.setTemperatureDelta(0.5);

    for (int i = 0; i < 1000; i++) {
        double oldTemp = t.calculateCurrentTemperature();
        t.apply();
        double newTemp = t.calculateCurrentTemperature();
        EXPECT_TRUE(oldTemp >= newTemp);
    }
}

// Holding Temperature
TEST(Thermostat, HoldingTemperature)
{
    BoundaryCondition o = BoundaryCondition::Outflow;
    BoundaryCondition r = BoundaryCondition::Reflecting;
    std::array<double, 3> domain = {80.0, 80.0, 1.0};
    std::array<BoundaryCondition, 6> bound = {r, r, r, r, o, o};
    double cutoff = 3;
    int parallelization = 0; //since tests are compiled without OpenMP this would automatically reset to 0 anyways
    auto pc = std::make_shared<LinkedCellParticleContainer>(cutoff, domain, bound, parallelization);

    std::array<double, 3> x1 = {1.0, 1.0, 0.0};
    std::array<double, 3> x2 = {79.0, 79.0, 79.0};
    std::array<double, 3> v = {0.0, 0.0, 0.0};
    double m = 1;
    double sigma = 1;
    double epsilon = 5;
    
    pc->addParticle(x1, v, m, epsilon, sigma);
    pc->addParticle(x2, v, m, epsilon, sigma);

    Thermostat t = Thermostat(pc, 100.0, 2);
    t.setTargetTemperature(100.0);
    t.setTemperatureDelta(0.5);
    t.initializeBrownianMotion();

    for (int i = 0; i < 1000; i++) {
        double oldTemp = t.calculateCurrentTemperature();
        t.apply();
        double newTemp = t.calculateCurrentTemperature();
        EXPECT_TRUE(std::abs(newTemp - oldTemp) < 1.0);
    }
}