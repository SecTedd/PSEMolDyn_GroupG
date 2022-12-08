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
    domain = {180, 90, 1};
    std::array<double, 3> dummyDomain = {180,90,1};
    BoundaryCondition b = BoundaryCondition::Outflow;
    boundaries = {b,b,b,b,b,b};
    particleContainer.reset(new LinkedCellParticleContainer(std::pow(2, 1.0/6), 3.0, dummyDomain, boundaries));
    end_time = 100;
    delta_t = 0.014;
    sigma = 1;
    epsilon = 5;
    cutoff = 3;
    writeFrequency = 10; 
    baseName = "outputVTK";
    benchmark_iterations = 0;
    showMenu = false;
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
const void ProgramParameters::setSigma(double sigma) { this->sigma = sigma; }
const void ProgramParameters::setEpsilon(double epsilon) { this->epsilon = epsilon; }
const void ProgramParameters::setCutoff(double cutoff) { this->cutoff = cutoff; }
const void ProgramParameters::setDomain(std::array<int, 3> domain) { this->domain = domain; }
const void ProgramParameters::setBoundaries(std::array<BoundaryCondition, 6> boundaries) { this->boundaries = boundaries; }
const void ProgramParameters::setWriteFrequency(int writeFrequency) { this->writeFrequency = writeFrequency; }
const void ProgramParameters::setBaseName(std::string baseName) {this->baseName = baseName; }
const void ProgramParameters::setShowMenu(bool show_menu) { this->showMenu = show_menu; }
const int ProgramParameters::getBenchmarkIterations() const { return benchmark_iterations; }
std::shared_ptr<ParticleContainer> ProgramParameters::getParticleContainer() { return particleContainer; }
const double ProgramParameters::getEndTime() const { return end_time; }
const double ProgramParameters::getDeltaT() const { return delta_t; }
const double ProgramParameters::getSigma() const { return sigma; }
const double ProgramParameters::getEpsilon() const { return epsilon; }
const double ProgramParameters::getCutoff() const { return cutoff; }
const std::array<int, 3> ProgramParameters::getDomain() const { return domain; }
const std::array<BoundaryCondition, 6> ProgramParameters::getBoundaries() const { return boundaries; }
const int ProgramParameters::getWriteFrequency() { return writeFrequency; }
const std::string ProgramParameters::getBaseName() { return baseName; }
const bool ProgramParameters::getShowMenu() const { return showMenu; }
