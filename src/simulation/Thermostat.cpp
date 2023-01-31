/*
 *  Thermostat.cpp
 *
 *  Created on: 14.12.2022
 *      Author: borisov
 */

#include "./Thermostat.h"
#include "../utils/ArrayUtils.h"
#include "../utils/MaxwellBoltzmannDistribution.h"


// Constructors

Thermostat::Thermostat(std::shared_ptr<ParticleContainer> particleContainer, double initTemperature, int dimension) {
    // Spdlog
    _logicLogger = spdlog::get("simulation_logger");
    _memoryLogger = spdlog::get("memory_logger");

    this->particleContainer = particleContainer;
    if (initTemperature < 0) {
        _logicLogger->info("Thermostat: Init temperature must be positive or zero!");
        _logicLogger->info("Thermostat: Using the absolute value of the provided init temperature!");
        initTemperature *= -1;
    }
    this->initTemperature = initTemperature;
    this->targetTemperature = initTemperature;
    this->temperatureDelta = -1;
    this->dimension = dimension;
    // default values: thermostat is applied to all 3 dimensions and all 3 dimensions contribute to temperature
    this->applyTo = {1, 1, 1};

    _memoryLogger->info("Thermostat generated!");
}

// Destructor

Thermostat::~Thermostat() {
    _memoryLogger->info("Thermostat destructed!");
}

void Thermostat::apply() {
    double currentTemperature = calculateCurrentTemperature();
    double newTemperature = calculateNewTemperature(currentTemperature);

    // skip when is exactly zero, to not divide by zero later
    if (currentTemperature == 0) {
        return;
    }

    // calculate beta
    double beta = sqrt(newTemperature / currentTemperature);

    for (auto &p: particleContainer->getActiveParticles()) {
        std::array<double, 3> velocity = p.getV();
        for (int i = 0; i < dimension; i++) {
            if (applyTo[i] == 1) {
                velocity[i] *= beta;
            }
        }
        p.setV(velocity);
    }

    _logicLogger->info("Temperature set to {} on directions {} {} {}", newTemperature, applyTo[0], applyTo[1], applyTo[2]);
}

void Thermostat::initializeBrownianMotion() {

    for (auto &p: particleContainer->getActiveParticles()) {
        double meanV = sqrt(this->initTemperature / p.getM());
        p.setV(p.getV() + maxwellBoltzmannDistributedVelocity(meanV, dimension));
    }

    _logicLogger->info("Initialized brownian motion with initial temperature: {}", this->initTemperature);
}

double Thermostat::calculateCurrentTemperature() {
    // calculate the kinetic energie
    double kineticE = 0;

    if (particleContainer->size() == 0) {
        return 0;
    }

    // subtract directions wich should not contribute to temperature
    std::array<double, 3> meanVs = {0, 0, 0};

    for (auto &p: particleContainer->getActiveParticles()) {
        meanVs[0] += p.getV()[0];
        meanVs[1] += p.getV()[1];
        meanVs[2] += p.getV()[2];
    }
    meanVs[0] /= particleContainer->size();
    meanVs[1] /= particleContainer->size();
    meanVs[2] /= particleContainer->size();

    for (auto &p: particleContainer->getActiveParticles()) {
        double dotProduct = 0;
        std::array<double, 3> velocity = p.getV();
        for (int j = 0; j < dimension; j++) {
            if (applyTo[j] == 0) {
                velocity[j] -= meanVs[j];
            }
        }
        for (long unsigned int i = 0; i < velocity.size(); i++) {
            dotProduct += velocity[i] * velocity[i];
        }
        kineticE += (p.getM() * dotProduct) / 2;
    }

    // calculate temperature from kinetic energie
    double temperature = (2 * kineticE) / (dimension * particleContainer->size());

    _logicLogger->debug("Current temperature: {}", temperature);

    return temperature;
}

double Thermostat::calculateNewTemperature(double currentTemperature) {
    // no deltaTemperature set
    if (temperatureDelta == -1) {
        return targetTemperature;
    }

    // deltaTemperature set

    // deltaTemperature is bigger than what needs to be changed
    if (std::abs(targetTemperature - currentTemperature) <= temperatureDelta) {
        return targetTemperature;
    }

    // current temperature is too high
    if (currentTemperature > targetTemperature) {
        return currentTemperature - temperatureDelta;
    }
    // current temperature is too low
    if (currentTemperature < targetTemperature) {
        return currentTemperature + temperatureDelta;
    }
    // curr == targetTemperature
    return currentTemperature;
}

// Getters

const double Thermostat::getTargetTemperature() {
    return this->targetTemperature;
}

const double Thermostat::getTemperatureDelta() {
    return this->temperatureDelta;
}

const double Thermostat::getInitTemperature() {
    return this->initTemperature;
}

const std::array<int, 3> Thermostat::getApplyTo() {
    return this->applyTo;
}

// Setters

const void Thermostat::setTargetTemperature(double targetTemperature) {
    if (targetTemperature < 0) {
        _logicLogger->info("Thermostat: Target temperature must be positive or zero!");
        _logicLogger->info("Thermostat: Using the absolute value of the provided target temperature!");
        targetTemperature *= -1;
    }
    this->targetTemperature = targetTemperature;
}

const void Thermostat::setTemperatureDelta(double temperatureDelta) {
    if (temperatureDelta < 0) {
        _logicLogger->info("Thermostat: Delta temperature must be positive or zero!");
        _logicLogger->info("Thermostat: Using the absolute value of the provided delta temperature!");
        temperatureDelta *= -1;
    }
    this->temperatureDelta = temperatureDelta;
}

const void Thermostat::setApplyTo(std::array<int, 3> applyTo) {
    this->applyTo = applyTo;
}