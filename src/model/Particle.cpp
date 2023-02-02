/*
 * Particle.cpp
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#include "Particle.h"
#include "../utils/ArrayUtils.h"

#include <iostream>
#include <omp.h>

Particle::Particle(int type_arg)
{
    type = type_arg;
    f = {0., 0., 0.};
    old_f = {0., 0., 0.};
    cell_idx = 0;
    invalid = false;
    halo = false;
    fixed = false;
    stiffness = 1;
    averageBondLength = 1;
    _memoryLogger = spdlog::get("memory_logger");
    _memoryLogger->info("Particle generated!");
}

Particle::Particle(const Particle &other)
{
    x = other.x;
    v = other.v;
    f = other.f;
    old_f = other.old_f;
    m = other.m;
    cell_idx = other.cell_idx;
    invalid = other.invalid;
    halo = other.halo;
    epsilon = other.epsilon;
    fixed = other.fixed;
    sigma = other.sigma;
    type = other.type;
    stiffness = other.stiffness;
    averageBondLength = other.averageBondLength;
    parallelNeighbours = other.parallelNeighbours; 
    diagonalNeighbours = other.diagonalNeighbours; 
    _memoryLogger = spdlog::get("memory_logger");
    _memoryLogger->info("Particle generated by copy!");
}

Particle::Particle(std::array<double, 3> x_arg, std::array<double, 3> v_arg,
                   double m_arg, double epsilon_arg, double sigma_arg, int type_arg, double stiffness_arg, double averageBondLength_arg) : x(x_arg), v(v_arg), m(m_arg), epsilon(epsilon_arg), sigma(sigma_arg), type(type_arg), stiffness(stiffness_arg), averageBondLength(averageBondLength_arg)
{
    f = {0., 0., 0.};
    old_f = {0., 0., 0.};
    cell_idx = 0;
    invalid = false;
    halo = false;
    fixed = false;
    _memoryLogger = spdlog::get("memory_logger");
    _memoryLogger->info("Particle generated!");
}

Particle::Particle(std::array<double, 3> x_arg, std::array<double, 3> v_arg,
                   double m_arg, double epsilon_arg, double sigma_arg, int type_arg, double stiffness_arg, double averageBondLength_arg, bool fixed_arg) : x(x_arg), v(v_arg), m(m_arg), epsilon(epsilon_arg), sigma(sigma_arg), type(type_arg), stiffness(stiffness_arg), averageBondLength(averageBondLength_arg), fixed(fixed_arg)
{
    f = {0., 0., 0.};
    old_f = {0., 0., 0.};
    cell_idx = 0;
    invalid = false;
    halo = false;
    _memoryLogger = spdlog::get("memory_logger");
    _memoryLogger->info("Particle generated!");
}

Particle::Particle(std::array<double, 3> x_arg, std::array<double, 3> v_arg,
                   double m_arg, double epsilon_arg, double sigma_arg, int type_arg) : x(x_arg), v(v_arg), m(m_arg), epsilon(epsilon_arg), sigma(sigma_arg), type(type_arg)
{
    f = {0., 0., 0.};
    old_f = {0., 0., 0.};
    cell_idx = 0;
    invalid = false;
    halo = false;
    fixed = false;
    stiffness = 1;
    averageBondLength = 1;
    _memoryLogger = spdlog::get("memory_logger");
    _memoryLogger->info("Particle generated!");
}

Particle::Particle(std::array<double, 3> x_arg, std::array<double, 3> v_arg,
                   double m_arg, double epsilon_arg, double sigma_arg, bool fixed_arg, int type_arg) : x(x_arg), v(v_arg), m(m_arg), epsilon(epsilon_arg), sigma(sigma_arg), type(type_arg), fixed(fixed_arg)
{
    f = {0.0, 0.0, 0.0};
    old_f = {0.0, 0.0, 0.0};
    cell_idx = 0;
    invalid = false;
    halo = false;
    stiffness = 1;
    averageBondLength = 1;
    _memoryLogger = spdlog::get("memory_logger");
    _memoryLogger->info("Particle generated!");
}

Particle::Particle(std::array<double, 3> x_arg, std::array<double, 3> v_arg, std::array<double, 3> f_arg, std::array<double, 3> old_f_arg,
                   double m_arg, double epsilon_arg, double sigma_arg, int type_arg, double stiffness_arg, double averageBondLength_arg) : x(x_arg), v(v_arg), f(f_arg), old_f(old_f_arg), m(m_arg), epsilon(epsilon_arg), sigma(sigma_arg), type(type_arg), stiffness(stiffness_arg), averageBondLength(averageBondLength_arg)
{
    cell_idx = 0;
    invalid = false;
    halo = false;
    fixed = false;
    _memoryLogger = spdlog::get("memory_logger");
    _memoryLogger->info("Particle generated!");
}

Particle::~Particle()
{
    _memoryLogger->info("Particle destructed!");
}

const std::array<double, 3> &Particle::getX() const { return x; }
const void Particle::setX(const std::array<double, 3> &new_x) {
    if (fixed) {
        return; 
    }
    x = new_x; 
}

const std::array<double, 3> &Particle::getV() const { return v; }
const void Particle::setV(const std::array<double, 3> &new_v) {
    if (fixed) {
        return;
    }
    v = new_v;
}

const std::array<double, 3> &Particle::getF() const { return f; }
const void Particle::setF(const std::array<double, 3> &new_f) {
    if (fixed) {
        return;
    }
    f = new_f;
}
const void Particle::addF(const std::array<double, 3> &new_f) {
    if (fixed) {
        return;
    }
    f = f + new_f;
}

const std::array<double, 3> &Particle::getOldF() const { return old_f; }
const void Particle::setOldF(const std::array<double, 3> &new_old_f) {
    if (fixed) {
        return;
    }
    old_f = new_old_f;
}

double Particle::getM() const { return m; }

int Particle::getType() const { return type; }

const int Particle::getCellIdx() { return cell_idx; }
const void Particle::setCellIdx(int cell_idx_arg) { cell_idx = cell_idx_arg; }

const bool Particle::getInvalid() { return invalid; }
const void Particle::setInvalid(bool invalid_arg) { invalid = invalid_arg; }

const bool Particle::getHalo() { return halo; }
const void Particle::setHalo(bool halo_arg) { halo = halo_arg; }

const double Particle::getEpsilon() const { return epsilon; }
const void Particle::setEpsilon(double epsilon_arg) { epsilon = epsilon_arg; }

const double Particle::getSigma() const { return sigma; }
const void Particle::setSigma(double sigma_arg) { sigma = sigma_arg; }

const bool Particle::getFixed() { return fixed; }
const void Particle::setFixed(bool fixed_arg) { fixed = fixed_arg; }

const double Particle::getStiffness() { return this->stiffness; }
const void Particle::setStiffness(double stiffness_arg) { this->stiffness = stiffness_arg; }

const double Particle::getAverageBondLength() { return this->averageBondLength; }
const void Particle::setAverageBondLength(double averageBondLength_arg) { this->averageBondLength = averageBondLength_arg; }

const std::vector<int> Particle::getParallelNeighbours() { return this->parallelNeighbours; }
const void Particle::setParallelNeighbours(std::vector<int> parallelNeighbours) { this->parallelNeighbours = parallelNeighbours; }

const std::vector<int> Particle::getDiagonalNeighbours() { return this->diagonalNeighbours; }
const void Particle::setDiagonalNeighbours(std::vector<int> diagonalNeighbours) { this->diagonalNeighbours = diagonalNeighbours; }


std::string Particle::toString() const
{
    std::stringstream stream;
    stream << "Particle: X:" << x << " v: " << v << " f: " << f
           << " old_f: " << old_f << " type: " << type << " cell index: " << cell_idx << " invalid: " << invalid << " halo: " << halo;
    return stream.str();
}

bool Particle::operator==(Particle &other)
{
    return (x == other.x) and (v == other.v) and (f == other.f) and
           (type == other.type) and (m == other.m) and (old_f == other.old_f);
}

std::ostream &operator<<(std::ostream &stream, Particle &p)
{
    stream << p.toString();
    return stream;
}
