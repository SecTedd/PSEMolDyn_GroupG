/*
 * ProgramParameters.cpp
 *
 *  Created on: 17.11.2022
 *      Author: wohlrapp
 */

#include "ProgramParameters.h"
#include "DirectSumParticleContainer.h"
#include "LinkedCellParticleContainer.h"
#include "ParticleCell.h"
#include "spdlog/spdlog.h"

#include <iostream>
#include <cmath>

ProgramParameters::ProgramParameters()
{
    domain = {3, 3, 1};
    BoundaryCondition o = BoundaryCondition::Outflow;
    boundaries = {o, o, o, o, o, o};
    dimension = 2;
    end_time = 1;
    delta_t = 0.0005;
    cutoff = 3;
    parallel = 1;
    writeFrequency = 50;
    particleContainer.reset(new LinkedCellParticleContainer(cutoff, domain, boundaries, parallel));
    baseName = "outputVTK";
    temp_init = 40;
    brownianMotion = true;
    n_thermostats = 1000;
    // has to be -1!
    temp_target = -1;
    // has to be -1!
    delta_temp = -1;
    g_grav = std::array<double, 3>{0.0, -12.44, 0.0};
    benchmark_iterations = 0;
    showMenu = false;
    createCheckpoint = false; 
    thermostat_applyTo = {1, 1, 1};
    csv_writeFrequency = 0;
    num_bins = 50;
    membrane = false;
    forces = std::list<std::shared_ptr<SingleParticleForce>>();
    memoryLogger = spdlog::get("memory_logger");
    memoryLogger->info("ProgramParameters generated!");
}

ProgramParameters::~ProgramParameters()
{
    memoryLogger->info("ProgramParameters destructed!");
}

const void ProgramParameters::resetParameters()
{
    particleContainer->resetParticles();
}

const void ProgramParameters::setEndTime(double end_time) { this->end_time = end_time; }
const void ProgramParameters::setDeltaT(double delta_t) { this->delta_t = delta_t; }
const void ProgramParameters::setBenchmarkIterations(int iterations) { this->benchmark_iterations = iterations; }
const void ProgramParameters::setCutoff(double cutoff)
{
    this->cutoff = cutoff;
    if (typeid(*particleContainer) == typeid(LinkedCellParticleContainer))
    {
        particleContainer.reset(new LinkedCellParticleContainer(cutoff, domain, boundaries, parallel));
    }
}
const void ProgramParameters::setDomain(std::array<double, 3> domain)
{
    this->domain = domain;
    if (typeid(*particleContainer) == typeid(LinkedCellParticleContainer))
    {
        particleContainer.reset(new LinkedCellParticleContainer(cutoff, domain, boundaries, parallel));
    }
}
const void ProgramParameters::setBoundaries(std::array<BoundaryCondition, 6> boundaries)
{
    this->boundaries = boundaries;

    // if 2D overwrite z-boundaries to be outflow
    if (this->domain[2] == 1)
    {
        this->boundaries[4] = BoundaryCondition::Outflow;
        this->boundaries[5] = BoundaryCondition::Outflow;
    }
    // error if periodic bondaries are not on opposite sides
    BoundaryCondition p = BoundaryCondition::Periodic;
    if ((boundaries[0] == p) != (boundaries[1] == p) || (boundaries[2] == p) != (boundaries[3] == p) || (boundaries[4] == p) != (boundaries[5] == p))
    {
        throw std::invalid_argument("Periodic boundaries have to be on opposite sides");
    }

    if (typeid(*particleContainer) == typeid(LinkedCellParticleContainer))
    {
        particleContainer.reset(new LinkedCellParticleContainer(cutoff, domain, boundaries, parallel));
    }
}

const void ProgramParameters::setParallel(int parallel) { 
    this->parallel = parallel; 
    if (typeid(*particleContainer) == typeid(LinkedCellParticleContainer))
     {
         particleContainer.reset(new LinkedCellParticleContainer(cutoff, domain, boundaries, parallel));
     }
}

const void ProgramParameters::setDimension(int dimension) { this->dimension = dimension; }
const void ProgramParameters::setWriteFrequency(int writeFrequency) { this->writeFrequency = writeFrequency; }
const void ProgramParameters::setBaseName(std::string baseName) { this->baseName = baseName; }
const void ProgramParameters::setTempInit(double temp_init) { this->temp_init = temp_init; }
const void ProgramParameters::setBrownianMotion(bool brownianMotion) { this->brownianMotion = brownianMotion; }
const void ProgramParameters::setNThermostats(int n_thermostats) { this->n_thermostats = n_thermostats; }
const void ProgramParameters::setTempTarget(double temp_target) { this->temp_target = temp_target; }
const void ProgramParameters::setDeltaTemp(double delta_temp) { this->delta_temp = delta_temp; }
const void ProgramParameters::setGGrav(std::array<double, 3> g_grav) { this->g_grav = g_grav; }
const void ProgramParameters::setShowMenu(bool show_menu) { this->showMenu = show_menu; }
const void ProgramParameters::setCreateCheckpoint(bool createCheckpoint) { this->createCheckpoint = createCheckpoint; }
const void ProgramParameters::setThermostatApplyTo(std::array<int, 3> thermostat_applyTo) { this->thermostat_applyTo = thermostat_applyTo; }
const void ProgramParameters::setCsvWriteFrequency(int csv_writeFrequency) { this->csv_writeFrequency = csv_writeFrequency; }
const void ProgramParameters::setNumBins(int num_bins) { this->num_bins = num_bins; }
const void ProgramParameters::setMembrane(bool membrane) { this->membrane = membrane; }
const void ProgramParameters::addForce(std::shared_ptr<SingleParticleForce> force) { forces.emplace_back(force); }
const int ProgramParameters::getBenchmarkIterations() const { return benchmark_iterations; }
std::shared_ptr<ParticleContainer> ProgramParameters::getParticleContainer() { return particleContainer; }
const double ProgramParameters::getEndTime() const { return end_time; }
const double ProgramParameters::getDeltaT() const { return delta_t; }
const double ProgramParameters::getCutoff() const { return cutoff; }
const std::array<double, 3> ProgramParameters::getDomain() const { return domain; }
const int ProgramParameters::getDimension() const { return dimension; }
const std::array<BoundaryCondition, 6> ProgramParameters::getBoundaries() const { return boundaries; }
const int ProgramParameters::getWriteFrequency() { return writeFrequency; }
const double ProgramParameters::getTempInit() const { return temp_init; }
const bool ProgramParameters::getBrownianMotion() const { return brownianMotion; }
const int ProgramParameters::getNThermostats() const { return n_thermostats; }
const double ProgramParameters::getTempTarget() const { return temp_target; }
const double ProgramParameters::getDeltaTemp() const { return delta_temp; }
const std::array<double, 3> ProgramParameters::getGGrav() const { return g_grav; }
const std::string ProgramParameters::getBaseName() { return baseName; }
const bool ProgramParameters::getShowMenu() const { return showMenu; }
const bool ProgramParameters::getCreateCheckpoint() { return createCheckpoint; }
const int ProgramParameters::getParallel() { return parallel; }
const std::array<int, 3> ProgramParameters::getThermostatApplyTo() const { return thermostat_applyTo; }
const int ProgramParameters::getCsvWriteFrequency() const { return this->csv_writeFrequency; }
const int ProgramParameters::getNumBins() const { return this->num_bins; }
const bool ProgramParameters::getMembrane() { return membrane; }
const std::list<std::shared_ptr<SingleParticleForce>> ProgramParameters::getForces() { return forces; }
