/*
 * ProgramParameters.h
 *
 *  Created on: 17.11.2022
 *      Author: wohlrapp
 */

#pragma once

#include "../model/ParticleContainer.h"
#include "../model/DirectSumParticleContainer.h"
#include "../model/ParticleCell.h"
#include "../simulation/SingleParticleForce.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/spdlog.h"

#include <list>

class InputFacade;

/**
 * @brief wrapper for all program parameters. Makes it easy to safe and rerun with the same parameters
 */
class ProgramParameters
{
private:
    std::shared_ptr<ParticleContainer> particleContainer; /// container for all the particles
    double end_time;                                      /// end_time of the simulation
    double delta_t;                                       /// increase in step size for the time
    double cutoff;                                        /// cutoff for the linked cell algorith
    std::array<double, 3> domain;                         /// the size of the domain
    int dimension;                                       /// the number of dimensions (2 or 3)
    std::array<BoundaryCondition, 6> boundaries;          /// the boundaries for the simulation
    int writeFrequency;                                   /// the number of iterations after which an vtk file is written
    std::string baseName;                                 /// the path to the output folder
    double temp_init;                                     /// the initial temperature
    bool brownianMotion;                                  /// specifies if particles should be initialized with brownian motion
    int n_thermostats;                                    /// the number of iterations after which the thermostat is applied
    double temp_target;                                   /// the target temperature of the simulation
    double delta_temp;                                    /// the maximum increase in the temperature per iteration
    std::array<double, 3> g_grav;                         /// the gravitational constant for the simulation
    int benchmark_iterations;                             /// number of runs in benchmark mode, 0 for normal simulations
    bool showMenu;                                        /// true if menu should be shown, false otherwise
    bool createCheckpoint;                                /// true if a checkpoint should be created, false otherwise
    int parallel;                                         /// 0 for no parallelization, 1 for first parallel strategy, 2 for the other
    std::shared_ptr<spdlog::logger> memoryLogger;         /// a speedlog logger which logs construction and destruction of particles
    std::array<int, 3> thermostat_applyTo;                /// indicator array on wich directions the thermostat will be applied
    int csv_writeFrequency;                               /// the number of iterations after wich a csv file is written (0 means never)
    int num_bins;                                         /// number of bins in x direction for density and velocity calculation
    bool membrane;                                          /// true if the simulation should calculate force according to a membrane
    std::list<std::shared_ptr<SingleParticleForce>> forces; /// list of forces that are applied

public:
    /**
     * @brief constructor for the ProgramParameters, initialises all of the parameters of the class
     */
    ProgramParameters();
    ~ProgramParameters();

    /**
     * @brief runs the simulation with the parameters that are currently set
     */
    const void runWithCurrentParameters();

    /**
     * @brief removes all particles from the simulation
     */
    const void resetParameters();

    /**
     * Getters/Setters
     */
    const void setEndTime(double end_time);

    const void setDeltaT(double delta_t);

    const void setBenchmarkIterations(int iterations);

    const void setCutoff(double cuttoff);

    const void setDomain(std::array<double, 3> domain);

    const void setDimension(int dimension);

    const void setBoundaries(std::array<BoundaryCondition, 6> boundaries);

    const void setWriteFrequency(int writeFrequency);

    const void setBaseName(std::string baseName);

    const void setTempInit(double temp_init);

    const void setBrownianMotion(bool browninanMotion);

    const void setNThermostats(int n_thermostats);

    const void setTempTarget(double temp_target);

    const void setDeltaTemp(double delta_temp);

    const void setGGrav(std::array<double, 3> g_grav);

    const void setShowMenu(bool show_menu);

    const void setCreateCheckpoint(bool createCheckpoint);

    const void setParallel(int parallel);
 
    const void setThermostatApplyTo(std::array<int, 3> thermostat_applyTo);

    const void setCsvWriteFrequency(int csv_writeFrequency);

    const void setNumBins(int num_bins);
    
    const void setMembrane(bool membrane);

    const void addForce(std::shared_ptr<SingleParticleForce> force); 

    std::shared_ptr<ParticleContainer> getParticleContainer();

    const double getEndTime() const;

    const double getDeltaT() const;

    const int getBenchmarkIterations() const;

    const double getCutoff() const;

    const std::array<double, 3> getDomain() const;

    const int getDimension() const;

    const std::array<BoundaryCondition, 6> getBoundaries() const;

    const int getWriteFrequency();

    const std::string getBaseName();

    const double getTempInit() const;

    const bool getBrownianMotion() const;

    const int getNThermostats() const;

    const double getTempTarget() const;

    const double getDeltaTemp() const;

    const std::array<double, 3> getGGrav() const;

    const bool getShowMenu() const;

    const bool getCreateCheckpoint(); 

    const int getParallel();
    
    const std::array<int, 3> getThermostatApplyTo() const;

    const int getCsvWriteFrequency() const;

    const int getNumBins() const;

    const bool getMembrane();

    const std::list<std::shared_ptr<SingleParticleForce>> getForces(); 
};