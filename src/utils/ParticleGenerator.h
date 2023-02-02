/*
 * ParticleGenerator.h
 *
 *  Created on: 04.12.2022
 *      Author: wohlrapp
 */

#pragma once

#include "../model/ParticleContainer.h"
#include "../model/Cuboid.h"
#include "../model/Sphere.h"
#include "./ArrayUtils.h"

#include <iostream>

namespace ParticleGenerator
{

    int index3DTo1D(std::array<int, 3> index, std::array<int, 3> n);
    std::array<int, 3> index1DTo3D(int index, std::array<int, 3> n);
    std::vector<int> getDiagonalNeighbourIndices(std::array<int, 3> n, double index);
    std::vector<int> getParallelNeighbourIndices(std::array<int, 3> n, double index);

    /**
     * @brief generates all particles of a cuboid and adds them to the particle container
     * @param particleContainer contains all particles for the simulation
     * @param cuboid contains all necessary parameters for the construction of the cuboid
     */
    inline void generateCuboid(ParticleContainer &particleContainer, Cuboid &cuboid, bool membrane)
    {
        // Variable init
        std::array<double, 3> lowerLeftCorner = cuboid.getX();
        std::array<double, 3> initV = cuboid.getV();
        double m = cuboid.getM();
        double epsilon = cuboid.getEpsilon();
        double sigma = cuboid.getSigma();
        int type = cuboid.getType();
        bool fixed = cuboid.getFixed();
        double stiffness = cuboid.getStiffness(); 
        double averageBondLength = cuboid.getAverageBondLength();
        std::array<int, 3> n = cuboid.getN();
        int numParticles = n[0] * n[1] * n[2];

        double meshWidth = cuboid.getH();

        // reserve memory for particles
        particleContainer.reserveMemoryForParticles(numParticles);

        // create particles
        std::array<double, 3> position;

        for (int x = 0; x < n[0]; x++)
        {
            for (int y = 0; y < n[1]; y++)
            {
                for (int z = 0; z < n[2]; z++)
                {
                    // initialize brownian motion
                    position[0] = lowerLeftCorner[0] + (x * meshWidth);
                    position[1] = lowerLeftCorner[1] + (y * meshWidth);
                    position[2] = lowerLeftCorner[2] + (z * meshWidth);
                    particleContainer.addParticle(position, initV, m, epsilon, sigma, type, stiffness, averageBondLength, fixed);
                }
            }
        }

        if (membrane)
        {
            auto &particles = particleContainer.getActiveParticles();
            int index = 0;
            for (int x = 0; x < n[0]; x++)
            {
                for (int y = 0; y < n[1]; y++)
                {
                    particles.at(index).setParallelNeighbours(getParallelNeighbourIndices(cuboid.getN(), index));
                    particles.at(index).setDiagonalNeighbours(getDiagonalNeighbourIndices(cuboid.getN(), index));
                    index++;
                }
            }
        }
    }

    inline std::vector<int> getParallelNeighbourIndices(std::array<int, 3> n, double index)
    {
        std::vector<std::array<int, 3>> neighbours3D;
        std::array<int, 3> index3D = index1DTo3D(index, n);

        int minX = index3D[0] - 1;
        int maxX = index3D[0] + 1;
        int minY = index3D[1] - 1;
        int maxY = index3D[1] + 1;

        if (index3D[0] == 0)
            minX = 0;
        if (index3D[0] == n[0] - 1)
            maxX = index3D[0];
        if (index3D[1] == 0)
            minY = 0;
        if (index3D[1] == n[1] - 1)
            maxY = index3D[1];

        for (int x = minX; x <= maxX; x++)
        {
            for (int y = minY; y <= maxY; y++)
            {
                if (y == index3D[1] || x == index3D[0])
                    neighbours3D.push_back(std::array<int, 3>{x, y, index3D[2]});
            }
        }

        std::vector<int> neighbours1D;

        for (auto i : neighbours3D)
        {
            int currentIndex = index3DTo1D(i, n);
            if (currentIndex > index)
                neighbours1D.emplace_back(currentIndex);
        }

        return neighbours1D;
    }

    inline std::vector<int> getDiagonalNeighbourIndices(std::array<int, 3> n, double index)
    {
        std::vector<std::array<int, 3>> neighbours3D;
        std::array<int, 3> index3D = index1DTo3D(index, n);

        int minX = index3D[0] - 1;
        int maxX = index3D[0] + 1;
        int minY = index3D[1] - 1;
        int maxY = index3D[1] + 1;

        if (index3D[0] == 0)
            minX = 0;
        if (index3D[0] == n[0] - 1)
            maxX = index3D[0];
        if (index3D[1] == 0)
            minY = 0;
        if (index3D[1] == n[1] - 1)
            maxY = index3D[1];

        for (int x = minX; x <= maxX; x++)
        {
            for (int y = minY; y <= maxY; y++)
            {
                if (y != index3D[1] && x != index3D[0])
                    neighbours3D.push_back(std::array<int, 3>{x, y, index3D[2]});
            }
        }

        std::vector<int> neighbours1D;

        for (auto i : neighbours3D)
        {
            int currentIndex = index3DTo1D(i, n);
            if (currentIndex > index)
                neighbours1D.emplace_back(currentIndex);
        }

        return neighbours1D;
    }

    inline std::array<int, 3> index1DTo3D(int index, std::array<int, 3> n)
    {
        std::array<int, 3> result;

        // x
        result[0] = index % n[0];
        // y
        result[1] = (index / n[0]) % n[1];
        // z
        result[2] = 1;

        return result;
    }

    inline int index3DTo1D(std::array<int, 3> index, std::array<int, 3> n)
    {
        return index[0] + (index[1] * n[0]);
    }

    /**
     * @brief generates all particles of a sphere and adds them to the particle container, if radius is 0, there is still
     * @param particleContainer contains all particles for the simulation
     * @param sphere contains all necessary parameters for the construction of the sphere
     * @param dimension 2D or 3D Sphere
     */
    inline void generateSphere(ParticleContainer &particleContainer, Sphere &sphere, int dimension)
    {
        std::array<double, 3> center = sphere.getCenter();
        double m = sphere.getM();
        int r = sphere.getR();
        double meshWidth = sphere.getH();
        std::array<double, 3> initV = sphere.getV();
        double epsilon = sphere.getEpsilon();
        double sigma = sphere.getSigma();
        int type = sphere.getType();
        bool fixed = sphere.getFixed();

        // number of particles which are later allocated
        int numParticles = 0;

        // create particles
        std::array<double, 3> position;

        // Variable init
        std::array<double, 3> startingPoint =
            {center[0] - (r - 1) * meshWidth - 0.5 * meshWidth, center[1] - (r - 1) * meshWidth - 0.5 * meshWidth, center[2] - ((dimension - 2) * ((r - 1) * meshWidth - 0.5 * meshWidth))};

        // first we need to iterate the loop once to count the number of particles
        for (int x = 0; x < 2 * r; x++)
        {
            for (int y = 0; y < 2 * r; y++)
            {
                position[0] = startingPoint[0] + (x * meshWidth);
                position[1] = startingPoint[1] + (y * meshWidth);
                position[2] = startingPoint[2];
                if (ArrayUtils::L2Norm(position - center) <= r * meshWidth)
                {
                    numParticles++;
                }

                // if three dimensions, we need to do this for every z
                for (int z = 1; z < (dimension - 2) * 2 * r; z++)
                {
                    position[0] = startingPoint[0] + (x * meshWidth);
                    position[1] = startingPoint[1] + (y * meshWidth);
                    position[2] = startingPoint[2] + (z * meshWidth);
                    if (ArrayUtils::L2Norm(position - center) <= r * meshWidth)
                    {
                        numParticles++;
                    }
                }
            }
        }

        // reserve memory for particles
        particleContainer.reserveMemoryForParticles(numParticles);

        // iterate again to actually add the particles to the container
        for (int x = 0; x < 2 * r; x++)
        {
            for (int y = 0; y < 2 * r; y++)
            {
                // do this once with z zero if we only have two dimensions
                // initialize brownian motion
                position[0] = startingPoint[0] + (x * meshWidth);
                position[1] = startingPoint[1] + (y * meshWidth);
                position[2] = startingPoint[2];

                // normally r-0.5 * mesh width but we want to include a bit more particles
                if (ArrayUtils::L2Norm(position - center) <= r * meshWidth)
                {
                    particleContainer.addParticle(position, initV, m, epsilon, sigma, fixed, type);
                }

                // if three dimensions, we need to do this for every z
                for (int z = 1; z < (dimension - 2) * 2 * r; z++)
                {

                    // initialize brownian motion
                    position[0] = startingPoint[0] + (x * meshWidth);
                    position[1] = startingPoint[1] + (y * meshWidth);
                    position[2] = startingPoint[2] + (z * meshWidth);
                    if (ArrayUtils::L2Norm(position - center) <= r * meshWidth)
                    {
                        particleContainer.addParticle(position, initV, m, epsilon, sigma, fixed, type);
                    }
                }
            }
        }
    }
}
