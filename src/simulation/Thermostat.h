/*
 *  Thermostat.h
 *
 *  Created on: 14.12.2022
 *      Author: borisov
 */

#pragma once

#include "../model/Particle.h"
#include "../model/ParticleContainer.h"
#include "spdlog/spdlog.h"

/**
 * @brief class that regulates the temperature within a system containing particles
 */
class Thermostat {
private:
    // desired temperature in wich the system wants to be
    double targetTemperature;

    // maximum absolute temperature change allowed for one application of the thermostat
    double temperatureDelta;

    // init temperature of the system
    double initTemperature;

    // dimensions of the observed system
    int dimension;

    // an array with either 1's or 0's indicating the directions where the velocity scaling will be applied to (e.g. only in x direction) 
    std::array<int, 3> applyTo;

    // an array indicating wich velocity direction will not contribute to the temperature (0 0 0 -> all directions contribute)
    std::array<int, 3> subtractMeanV;

    // particles of the observed system
    std::shared_ptr<ParticleContainer> particleContainer;
    
    // a speedlog logger which logs the logic flow of the simulation
    std::shared_ptr<spdlog::logger> _logicLogger;
    
    //a speedlog logger which logs construction and destruction of particles
    std::shared_ptr<spdlog::logger> _memoryLogger;

public:

    /**
     * @brief Construct a new Thermostat object
     * 
     * @param particleContainer particles of the observed system
     * @param initTemperature init temperature of the system the thermostat is watching
     * @param dimension dimensions of the observed system
     */
    Thermostat(std::shared_ptr<ParticleContainer> particleContainer, double initTemperature, int dimension);

    ~Thermostat();
    
    /**
     * @brief calculates the current Temperature of the system (using the kinetic energy)
     * 
     * @return double the current temperature of the system
     */
    double calculateCurrentTemperature();

    /**
     * @brief calculates the new temperature of a system, according to targetTemperature and temperatureDelta
     * 
     * @param currentTemperature current temperature in the system
     * @return double the new temperature the system should have after the thermostat is applied
     */
    double calculateNewTemperature(double currentTemperature);

    /**
     * @brief changes the temperature of the system towards the targetTemperature
     */
    void apply();

    /**
     * @brief initializes the velocity of the particles with the brownian motion
     */
    void initializeBrownianMotion();

    // Getters

    const double getTargetTemperature();

    const double getTemperatureDelta();

    const double getInitTemperature();

    const std::array<int, 3> getApplyTo();

    const std::array<int, 3> getSubtractMeanV();

    // Setters

    const void setTargetTemperature(double targetTemperature);

    const void setTemperatureDelta(double temperatureDelta);

    const void setApplyTo(std::array<int, 3> applyTo);

    const void setSubtractMeanV(std::array<int, 3> subtractMeanV);
};