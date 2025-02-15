/*
 * XYZWriter.h
 *
 *  Created on: 01.03.2010
 *      Author: eckhardw
 */

#pragma once

#include "../model/Particle.h"
#include "spdlog/spdlog.h"

#include <fstream>
#include <vector>
#include <memory>

namespace outputWriter
{

    class XYZWriter
    {
    private:
        /**
         * a speedlog logger which logs construction and destruction of particles
         */
        std::shared_ptr<spdlog::logger> _memoryLogger;

    public:
        XYZWriter();

        virtual ~XYZWriter();

        void plotParticles(std::vector<Particle> particles, const std::string &filename,
                           int iteration);
    };

} // namespace outputWriter
