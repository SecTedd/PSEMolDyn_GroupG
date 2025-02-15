/*
 * Simulation.cpp
 *
 * Created: 31.10.2022
 * Author:  wohlrapp
 */

#include "Simulation.h"
#include "../outputWriter/VTKWriter.h"
#include "../utils/ArrayUtils.h"
#include "../outputWriter/OutputFacade.h"
#include "LennardJonesForce.h"
#include "LennardJonesForceHarmonic.h"
#include "InterParticleGravitationalForce.h"
#include "SingleParticleGravitationalForce.h"
#include "../model/ProgramParameters.h"
#include "./Thermostat.h"
#include "./DVProfileCalculator.h"

#include <iostream>

Simulation::Simulation(ProgramParameters *programParameters)
{
    _programParameters = programParameters;
    if (programParameters->getMembrane())
        _interParticleForceCalculation.reset(new LennardJonesForceHarmonic());
    else
        _interParticleForceCalculation.reset(new LennardJonesForce());
    _singleParticleForceCalculations = programParameters->getForces();
    std::shared_ptr<SingleParticleForce> force; 
    force.reset(new SingleParticleGravitationalForce(programParameters->getGGrav()));
    _singleParticleForceCalculations.emplace_back(force);
    _logicLogger = spdlog::get("simulation_logger");
    _memoryLogger = spdlog::get("memory_logger");
    _memoryLogger->info("Simulation generated!");
}

Simulation::~Simulation()
{
    _memoryLogger->info("Simulation destructed!");
}

const void Simulation::simulate()
{
    double start_time = 0;
    double current_time = start_time;

    int iteration = 0;

    OutputFacade outputFacade = OutputFacade(_programParameters);

    // initialize Thermostat
    Thermostat t = Thermostat(_programParameters->getParticleContainer(), _programParameters->getTempInit(), _programParameters->getDimension());
    // target temperature provided
    if (_programParameters->getTempTarget() != -1)
    {
        t.setTargetTemperature(_programParameters->getTempTarget());
    }
    // temperature delta provided
    if (_programParameters->getDeltaTemp() != -1)
    {
        t.setTemperatureDelta(_programParameters->getDeltaTemp());
    }
    t.setApplyTo(_programParameters->getThermostatApplyTo());
    // initialize browninan motion if needed
    if (_programParameters->getBrownianMotion())
    {
        t.initializeBrownianMotion();
    }

    DVProfileCalculator dv_calc = DVProfileCalculator(_programParameters->getParticleContainer(), _programParameters->getNumBins(), _programParameters->getDomain());

    // calculating force once to initialize force
    _interParticleForceCalculation->calculateForce(*_programParameters->getParticleContainer());
    for (auto force : _singleParticleForceCalculations)
    {
        force->calculateForce(*_programParameters->getParticleContainer(), current_time);
    }

    //only write 0th file if not in benchmark mode
    if (_programParameters->getBenchmarkIterations() == 0)
    {
        outputFacade.outputVTK(iteration);
    }
        
    // for this loop, we assume: current x, current f and current v are known
    while (current_time < _programParameters->getEndTime())
    {
        // calculate new x
        calculateX();

        // calculate new f
        _interParticleForceCalculation->calculateForce(*_programParameters->getParticleContainer());
        for (auto force : _singleParticleForceCalculations)
        {
            force->calculateForce(*_programParameters->getParticleContainer(), current_time);
        }
        // calculate new v
        calculateV();

        // if n_thermostats = 0 the thermostat is off
        if (_programParameters->getNThermostats() != 0 && iteration % _programParameters->getNThermostats() == 0)
        {
            t.apply();
        }

        if (iteration % _programParameters->getWriteFrequency() == 0 && _programParameters->getBenchmarkIterations() == 0)
        {
            outputFacade.outputVTK(iteration);
        }

        // if csv_writeFrequency = 0 no csv output will be written
        if (_programParameters->getCsvWriteFrequency() != 0 && iteration % _programParameters->getCsvWriteFrequency() == 0 && _programParameters->getBenchmarkIterations() == 0) 
        {
            std::vector<int> data = dv_calc.calculate();
            outputFacade.writeCSV(data, dv_calc.getAvg());
        }
        _logicLogger->info("Iteration {} finished.", iteration);
        iteration++;

        current_time += _programParameters->getDeltaT();
    }
    if (_programParameters->getCreateCheckpoint())
        outputFacade.createCheckpoint();

    _logicLogger->info("Finished Iterations. Terminating");
}

void Simulation::calculateX()
{
    std::shared_ptr<ParticleContainer> particleContainer = _programParameters->getParticleContainer();

    // creating lambda to calculate new position based on the Velocity-Störmer-Verlet algortihm
    std::function<void(Particle &)> f = [delta_t = _programParameters->getDeltaT(), logicLogger = _logicLogger](Particle &p1)
    {
        std::array<double, 3> x_new = p1.getX() + delta_t * p1.getV() + (delta_t * delta_t / (2 * p1.getM())) * p1.getF();
        if (p1.getF()[0] >= 10e9 || p1.getF()[0] <= -10e9)
        {
            logicLogger->debug("High force: " + p1.toString());
            logicLogger->debug("New X: " + std::to_string(x_new[0]) + ", " + std::to_string(x_new[1]) + ", " + std::to_string(x_new[2]));
        }
        p1.setX(x_new);
    };

    particleContainer->iterateParticles(f, true);
}

void Simulation::calculateV()
{
    std::shared_ptr<ParticleContainer> particleContainer = _programParameters->getParticleContainer();

    // creating lambda to calculate new speed based on the Velocity-Störmer-Verlet algortihm
    std::function<void(Particle &)> f = [delta_t = _programParameters->getDeltaT()](Particle &p1)
    {
        std::array<double, 3> v_new = p1.getV() + (delta_t / (2 * p1.getM())) * (p1.getOldF() + p1.getF());
        p1.setV(v_new);
    };

    particleContainer->iterateParticles(f, false);
}

const std::shared_ptr<spdlog::logger> Simulation::getLogicLogger() const { return _logicLogger; }
