/*
 * Cuboid.h
 *
 *  Created on: 15.11.2022
 *      Author: borisov
 */

#pragma once

#include "spdlog/spdlog.h"
#include "../simulation/SingleParticleForce.h"

#include <array>
#include <list>

/**
 * @brief Cuboid data class, which holds all the information about a cuboid
 */
class Cuboid
{
private:
    // coordinate of the lower left front-side corner
    std::array<double, 3> x;

    // number of particles per dimension
    std::array<int, 3> n;

    // distance of the particles (mesh width of the grid)
    double h;

    // mass of one particle
    double m;

    // initial velocity of the particles
    std::array<double, 3> v;

    // type of the particles
    int type;

    // Variable for the Lennard Jones force
    double epsilon;

    // Variable for the Lennard Jones force
    double sigma;

    /**
     * Stiffness of molecule
     */
    double stiffness;

    /**
     * Average Bond Length of molecule
     */
    double averageBondLength;

    /**
     * a speedlog logger which logs construction and destruction of particles
     */
    std::shared_ptr<spdlog::logger> _memoryLogger;

public:
    /**
     * @brief Construct a new Cuboid object
     *
     * @param x The coordinate of the lower left front-side corner
     * @param n The number of particles per dimension
     * @param h The distance of the particles (mesh width of the grid)
     * @param m The mass of one particle
     * @param v The initial velocity of the particles
     * @param sigma The sigma value of the particles
     * @param epsilon The epsilon value of the particles in the cuboid
     * @param type The type of the particles in the cuboid
     * @param stiffness The stiffness of the molecule
     * @param averageBondLength The average Bond Length of the molecule
     */
    Cuboid(std::array<double, 3> x, std::array<int, 3> n, double h, double m, std::array<double, 3> v, double epsilon, double sigma, int type, double stiffness_arg = 1, double averageBondLength_arg = 1);

    ~Cuboid();

    /*
     * Getters
     */

    const std::array<double, 3> getX();

    const std::array<int, 3> getN();

    const double getH();

    const double getM();

    const std::array<double, 3> getV();

    const double getEpsilon() const;

    const double getSigma() const;

    const int getType();

    const double getStiffness();

    const double getAverageBondLength();

    /*
     * Setters
     */

    const void setX(std::array<double, 3> &x);

    const void setN(std::array<int, 3> &n);

    const void setH(double h);

    const void setM(double m);

    const void setV(std::array<double, 3> &v);

    const void setType(int type);

    const void setStiffness(double stiffness_arg);

    const void setAverageBondLength(double averageBondLength_arg);
};
