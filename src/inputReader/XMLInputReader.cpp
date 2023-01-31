/*
 * XMLInputReader.cpp
 *
 *  Created on: 07.12.2022
 *      Author: wohlrapp
 */

#include "./XMLInputReader.h"
#include "./InputFacade.h"
#include "../model/Sphere.h"
#include "../simulation/TemporalSingleParticleForce.h"
#include "../utils/MaxwellBoltzmannDistribution.h"
#include "../utils/ArrayUtils.h"
#include "../utils/ParticleGenerator.h"
#include "../xsd/Simulation.cxx"
#include "../xsd/SimulationState.cxx"
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

BoundaryCondition getBoundaryCondition(std::string s)
{
    if (s == "Reflecting")
        return BoundaryCondition::Reflecting;

    if (s == "Periodic")
        return BoundaryCondition::Periodic;

    return BoundaryCondition::Outflow;
}

void XMLInputReader::readInput(ProgramParameters &programParameters, const char *filename)
{

    // not the best solution but not possible otherwise
    try
    {
        std::shared_ptr<simulation_state_t> state(simulation_state(filename));

        for (simulation_state_t::particle_const_iterator i(state->particle().begin()); i != state->particle().end(); i++)
        {
            std::array<double, 3> position;
            simulation_state_t::particle_traits::type::x_type x_xml = i->x();
            position[0] = x_xml.x1();
            position[1] = x_xml.y();
            position[2] = x_xml.z();

            std::array<double, 3> velocity;
            simulation_state_t::particle_traits::type::v_type v = i->v();
            velocity[0] = v.x();
            velocity[1] = v.y();
            velocity[2] = v.z();

            double m = i->mass();
            double epsilon = i->epsilon();
            double sigma = i->sigma();

            std::array<double, 3> f;
            simulation_state_t::particle_traits::type::f_type f_xml = i->f();
            f[0] = f_xml.x();
            f[1] = f_xml.y();
            f[2] = f_xml.z();

            std::array<double, 3> old_f;
            simulation_state_t::particle_traits::type::old_f_type old_f_xml = i->old_f();
            old_f[0] = old_f_xml.x();
            old_f[1] = old_f_xml.y();
            old_f[2] = old_f_xml.z();

            int type = i->type();
            double stiffness = i->stiffness(); 
            double averageBondLength = i->averageBondLength();

            programParameters.getParticleContainer()->addParticle(position, velocity, f, old_f, m, epsilon, sigma, type, stiffness, averageBondLength);
        }
        return;
    }
    catch (xml_schema::unexpected_element const &)
    {
        // trying next function
        getLogicLogger()->info("Not SimulationState");
    }
    catch (const xml_schema::exception &e)
    {
        getLogicLogger()->error(e.what());
    }

    try
    {
        std::shared_ptr<simulation_t> xml(simulation(filename));

        std::shared_ptr<InputFacade> inputFacade = std::make_shared<InputFacade>();

        programParameters.setEndTime(xml->end_time());
        programParameters.setDeltaT(xml->delta_t());
        programParameters.setCutoff(xml->cutoff());

        if (xml->brownianMotion().present())
            programParameters.setBrownianMotion(xml->brownianMotion().get());

        std::array<double, 3> domain;
        simulation_t::domain_type d = xml->domain();
        domain[0] = d.x();
        domain[1] = d.y();
        domain[2] = d.z();
        programParameters.setDomain(domain);

        if (xml->dimension().present())
            programParameters.setDimension(xml->dimension().get());

        std::array<BoundaryCondition, 6> boundaries = std::array<BoundaryCondition, 6>();
        simulation_t::boundaries_type b = xml->boundaries();
        std::string boundary = b.xLeft();
        boundaries[0] = getBoundaryCondition(boundary);
        boundary = b.xRight();
        boundaries[1] = getBoundaryCondition(boundary);
        boundary = b.yBottom();
        boundaries[2] = getBoundaryCondition(boundary);
        boundary = b.yTop();
        boundaries[3] = getBoundaryCondition(boundary);
        boundary = b.zFront();
        boundaries[4] = getBoundaryCondition(boundary);
        boundary = b.zBack();
        boundaries[5] = getBoundaryCondition(boundary);
        programParameters.setBoundaries(boundaries);

        std::array<double, 3> gGrav;
        simulation_t::g_grav_type g = xml->g_grav();
        gGrav[0] = g.x();
        gGrav[1] = g.y();
        gGrav[2] = g.z();
        programParameters.setGGrav(gGrav);

        if (xml->membrane().present())
        {
            programParameters.setMembrane(xml->membrane().get());
        }

        if (xml->writeFrequency().present())
        {
            programParameters.setWriteFrequency(xml->writeFrequency().get());
        }

        if (xml->baseName().present())
        {
            programParameters.setBaseName(xml->baseName().get());
        }

        if (xml->createCheckpoint().present())
        {
            programParameters.setCreateCheckpoint(xml->createCheckpoint().get());
        }

        for (simulation_t::file_name_const_iterator i(xml->file_name().begin()); i != xml->file_name().end(); i++)
        {
            std::string filename = i->substr(0, i->length());
            inputFacade->readInput(programParameters, filename.c_str());
        }

        if(xml->csvWriteFrequency().present()){
            programParameters.setCsvWriteFrequency(xml->csvWriteFrequency().get()); 
        }

        if(xml->numBins().present()){
            programParameters.setNumBins(xml->numBins().get()); 
        }

        if (xml->thermostat().present())
        {
            auto thermostat = xml->thermostat().get();

            programParameters.setNThermostats(thermostat.n_thermostat());

            std::array<int, 3> applyTo;

            applyTo[0] = thermostat.apply_to().x();
            applyTo[1] = thermostat.apply_to().y();
            applyTo[2] = thermostat.apply_to().z();

            programParameters.setThermostatApplyTo(applyTo);

            if (thermostat.temp_init().present())
                programParameters.setTempInit(thermostat.temp_init().get());

            if (thermostat.temp_target().present())
            {
                programParameters.setTempTarget(thermostat.temp_target().get());
            }

            if (thermostat.delta_temp().present())
            {
                programParameters.setDeltaTemp(thermostat.delta_temp().get());
            }
        }

        for (simulation_t::cuboid_const_iterator i(xml->cuboid().begin()); i != xml->cuboid().end(); i++)
        {
            std::array<double, 3> position;
            simulation_t::cuboid_traits::type::position_type p = i->position();
            position[0] = p.x();
            position[1] = p.y();
            position[2] = p.z();

            std::array<double, 3> velocity;

            simulation_t::cuboid_traits::type::velocity_type v = i->velocity();

            velocity[0] = v.x();
            velocity[1] = v.y();
            velocity[2] = v.z();

            std::array<int, 3> dimensions;

            simulation_t::cuboid_traits::type::dimensions_type dim = i->dimensions();
            dimensions[0] = dim.x();
            dimensions[1] = dim.y();
            dimensions[2] = dim.z();

            for (simulation_t::cuboid_type::force_const_iterator j(i->force().begin()); j != i->force().end(); j++)
            {
                std::array<double, 3> force;
                simulation_t::cuboid_type::force_type::force1_type f1 = j->force1();
                force[0] = f1.x();
                force[1] = f1.y();
                force[2] = f1.z();

                double end_time = j->end_time();

                std::vector<int> indices;
                int numberOfParticles = programParameters.getParticleContainer()->size();

                for (simulation_t::cuboid_type::force_type::particles_const_iterator k(j->particles().begin()); k != j->particles().end(); k++)
                {
                    std::array<int, 3> index3D;
                    simulation_t::cuboid_type::force_type::particles_type::particle_index_type pIndex = k->particle_index();
                    index3D[0] = pIndex.x();
                    index3D[1] = pIndex.y();
                    index3D[2] = pIndex.z();

                    indices.push_back(ParticleGenerator::index3DTo1D(index3D, dimensions) + numberOfParticles);
                }
                std::shared_ptr<SingleParticleForce> spForce;
                spForce.reset(new TemporalSingleParticleForce(force, end_time, indices));
                programParameters.addForce(spForce);
            }

            double h = i->h();
            double m = i->mass();
            double epsilon = i->epsilon();
            double sigma = i->sigma();
            int type = i->type();
            bool fixed = false;

            if (i->fixed().present())
            {
                fixed = i->fixed().get();
            }

            double stiffness = 1;
            double averageBondLength = 1;

            if (i->stiffness().present())
            {
                stiffness = i->stiffness().get();
            }

            if (i->average_bond_length().present())
            {
                averageBondLength = i->average_bond_length().get();
            }

            std::unique_ptr<Cuboid> cuboid = std::make_unique<Cuboid>(Cuboid(position, dimensions, h, m, velocity, epsilon, sigma, type, stiffness, averageBondLength, fixed));
            ParticleGenerator::generateCuboid(*programParameters.getParticleContainer(), *cuboid, programParameters.getMembrane());
        }

        for (simulation_t::sphere_const_iterator i(xml->sphere().begin()); i != xml->sphere().end(); i++)
        {
            std::array<double, 3> center;
            simulation_t::sphere_traits::type::center_type c = i->center();
            center[0] = c.x();
            center[1] = c.y();
            center[2] = c.z();

            std::array<double, 3> velocity;
            simulation_t::sphere_traits::type::velocity_type v = i->velocity();
            velocity[0] = v.x();
            velocity[1] = v.y();
            velocity[2] = v.z();

            int r = i->r();
            double h = i->h();
            double m = i->mass();
            double epsilon = i->epsilon();
            double sigma = i->sigma();
            int type = i->type();
            bool fixed = false;

            if (i->fixed().present())
            {
                fixed = i->fixed().get();
            }

            std::unique_ptr<Sphere> sphere = std::make_unique<Sphere>(Sphere(center, r, h, m, velocity, epsilon, sigma, type, fixed));
            ParticleGenerator::generateSphere(*programParameters.getParticleContainer(), *sphere, programParameters.getDimension());
        }
    }
    catch (const xml_schema::exception &e)
    {
        getLogicLogger()->error(e.what());
    }
}