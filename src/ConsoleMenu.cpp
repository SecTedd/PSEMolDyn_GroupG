/*
 *  ConsoleMenu.cpp
 *
 *  Created on: 17.11.2022
 *      Author: wohlrapp
 */

#include "spdlog/spdlog.h"
#include "ConsoleMenu.h"
#include "./simulation/Simulation.h"

#include <set>
#include <sstream>
#include <iostream>
#include <regex>
#include <fstream>
#include <chrono>

ConsoleMenu::ConsoleMenu(ProgramParameters *programParameters, InputFacade *inputFacade)
{
    _inputFacade = inputFacade;
    _programParameters = programParameters;
    _memoryLogger = spdlog::get("memory_logger");
    _memoryLogger->info("ConsoleMenu generated!");
}

ConsoleMenu::~ConsoleMenu()
{
    _memoryLogger->info("ConsoleMenu destructed!");
}

const void ConsoleMenu::openMenu()
{
    printHelpMenu();
    while (1)
    {
        std::string command;
        std::cout << "MolSim Group G > Enter command: ";
        std::getline(std::cin, command);
        command = Input::ftrim(command);
        if (!verifyCommand(command))
        {
            std::cout << "MolSim Group G > Invalid command: " << command << std::endl;
            printHelpMenu();
        }
        else
        {
            char c = command[1];
            std::string parameter;

            switch (c)
            {
            case 'f':
            {
                int num_particles_before = _programParameters->getParticleContainer()->size();
                parameter = Input::trim(Input::trim(command.substr(2)), "\"");
                _inputFacade->readInput(*_programParameters, parameter.c_str());
                int num_particles_after = _programParameters->getParticleContainer()->size();
                int num_added_particles = num_particles_after - num_particles_before;
                std::cout << "MolSim Group G > " << num_added_particles << " particles were added to the simulation" << std::endl;
            }
            break;
            case 'e':
                parameter = Input::trim(command.substr(2));
                _programParameters->setEndTime(std::__cxx11::stod(parameter));
                std::cout << "MolSim Group G > End time was set to " << _programParameters->getEndTime() << std::endl;
                break;
            case 'd':
                parameter = Input::trim(command.substr(2));
                _programParameters->setDeltaT(std::__cxx11::stod(parameter));
                std::cout << "MolSim Group G > Delta time was set to " << _programParameters->getDeltaT() << std::endl;
                break;
            case 'c': 
                _programParameters->setCreateCheckpoint(!_programParameters->getCreateCheckpoint());
                break; 
            case 'i':
            {
                int num_particles = _programParameters->getParticleContainer()->size();
                std::cout << "MolSim Group G > Number of particles in the simulation: " << num_particles << std::endl;
                if (num_particles > 0)
                {
                    std::cout << "MolSim Group G > You can reset the number of particles to 0 with the -x command" << std::endl;
                }
                std::cout << "MolSim Group G > Current end time: " << _programParameters->getEndTime() << std::endl;
                std::cout << "MolSim Group G > Current delta time: " << _programParameters->getDeltaT() << std::endl;
                std::cout << "MolSim Group G > Current createCheckpoint: " << _programParameters->getCreateCheckpoint() << std::endl;
                std::string mode = _programParameters->getBenchmarkIterations() == 0 ? "Simulation" : "Benchmark";
                std::cout << "MolSim Group G > Current mode: " << mode << std::endl;
            }
            break;
            case 'x':
                _programParameters->resetParameters();
                std::cout << "MolSim Group G > All particles were removed from the simulation" << std::endl;
                break;
            case 'r':
            {
                if (_programParameters->getBenchmarkIterations() == 0)
                {
                    std::cout << "MolSim Group G > Running simulation ... " << std::endl;
                    std::unique_ptr<Simulation> simulation = std::make_unique<Simulation>(Simulation(_programParameters));
                    simulation->simulate();
                    std::cout << "MolSim Group G > ... Finished" << std::endl;
                }
                else
                {
                    std::cout << "MolSim Group G > Running benchmark ..." << std::endl;

                    using namespace std::chrono;

                    time_point<high_resolution_clock> start_point, end_point;
                    auto total_time = microseconds(0).count();

                    for (int i = 0; i < _programParameters->getBenchmarkIterations(); i++)
                    {
                        std::unique_ptr<Simulation> simulation = std::make_unique<Simulation>(Simulation(_programParameters));

                        start_point = high_resolution_clock::now();
                        simulation->simulate();
                        end_point = high_resolution_clock::now();

                        auto start = time_point_cast<microseconds>(start_point).time_since_epoch().count();
                        auto end = time_point_cast<microseconds>(end_point).time_since_epoch().count();

                        total_time += (end - start);
                    }
                    auto mean_time = total_time / _programParameters->getBenchmarkIterations();
                    std::cout << "MolSim Group G > Mean duration over " << _programParameters->getBenchmarkIterations() << " run(s): " << mean_time / 1000000.0 << " seconds" << std::endl;
                    std::cout << "MolSim Group G > ... Finished" << std::endl;
                }
            }
            break;
            case 'h':
                printHelpMenu();
                break;
            case 'q':
                std::cout << "MolSim Group G > Quitting ..." << std::endl;
                return;
            }
            printf("===============================================================================================================================================================\n");
        }
    }
}

const bool ConsoleMenu::verifyCommand(std::string command) const
{
    std::set<char> commands = {'f', 'c', 'e', 'd', 'x', 'r', 'h', 'q', 'i'};
    if (command.length() < 2 || command[0] != '-' || commands.count(command[1]) == 0)
    {
        return false;
    }
    std::set<char> commandsWithParameters = {'f', 't', 'd'};
    bool parameterTest;
    if (commandsWithParameters.count(command[1]) != 0)
    {
        switch (command[1])
        {
        case 'f':
        {
            std::ifstream test(Input::trim(Input::trim(command.substr(2)), "\""));
            if (!test)
            {
                std::cout << "MolSim Group G > Error: File does not exist" << std::endl;
                return false;
            }
        }
        break;
        case 't':
        case 'd':
        case 'y':
        case 's':
            parameterTest = Input::isDouble(Input::trim(command.substr(2)));
            if (!parameterTest)
            {
                std::cout << "MolSim Group G > Error: Parameter is not a double" << std::endl;
                return false;
            }
            break;
        }
    }
    return true;
}

const void ConsoleMenu::printHelpMenu() const
{
    printf("===============================================================================================================================================================\n");
    printf("MolSim Group G > Welcome to the menu of the MolSim application. In here, you can set parameters of the application, read in multiple inputs, and finally run it\n");
    printf("MolSim Group G > -f <filename> .......... The path to an input file. If not specified no particles are generated\n");
    printf("MolSim Group G > -e <end_time> .......... The end time of the simulation. If not specified, 100 is used\n");
    printf("MolSim Group G > -d <delta_t> ........... The size of the time steps in the simulation. If not specified 0.014 is used\n");
    printf("MolSim Group G > -c ..................... If a createCheckpoint is 'false' it is set to 'true', if it is 'true' it is set to 'false'\n");
    printf("MolSim Group G > -i ..................... Displays the currently set values in the simulation\n");
    printf("MolSim Group G > -x ..................... Deletes all particles from the simulation\n");
    printf("MolSim Group G > -r ..................... Run the program with the currently set values\n");
    printf("MolSim Group G > -h ..................... Help for the menu\n");
    printf("MolSim Group G > -q ..................... Quit the program\n");
    printf("===============================================================================================================================================================\n");
}
