/*
 *  DVProfileCalculator.cpp
 *
 *  Created on: 25.01.2023
 *      Author: borisov
 */

#include "./DVProfileCalculator.h"
#include "../utils/ArrayUtils.h"

#include <iostream>

DVProfileCalculator::DVProfileCalculator(std::shared_ptr<ParticleContainer> particleContainer, int n, std::array<double, 3> domain) {
    // Spdlog
    _logicLogger = spdlog::get("simulation_logger");
    _memoryLogger = spdlog::get("memory_logger");

    this->particleContainer = particleContainer;
    this->n = n;
    this->domain = domain;

    for (int i = 0; i < n; i++) {
        avg.push_back(0);
    }

    _memoryLogger->info("DVProfielCalculator generated!");
}

DVProfileCalculator::~DVProfileCalculator() {
    _memoryLogger->info("DVProfileCalculator destructed!");
}

const std::vector<int> DVProfileCalculator::calculate() {
    // create array of n bins
    std::vector<int> bins;

    for (int i = 0; i < n; i++) {
        bins.push_back(0);
    }

    double bin_width = domain[0] / n; 

    for (auto &p: particleContainer->getActiveParticles()) {
        int bin = p.getX()[0] / bin_width;
        bins.at(bin)++;
    }

    this->avg = this->avg + bins;
    return bins;
}

// Getters
const int DVProfileCalculator::getN() {
    return this->n;
}

const std::vector<int> DVProfileCalculator::getAvg() {
    return this->avg;
}

// Setters
const void DVProfileCalculator::setN(int n) {
    this->n = n;
}

const void DVProfileCalculator::setAvg(std::vector<int> avg) {
    this->avg = avg;
}

