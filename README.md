# MolSim PSEMolDyn_GroupG

Code for the practical course *PSE: Molecular Dynamics* by group G.

Members: <br />
Alexandra Marquardt <br />
Matteo Wohlrapp <br />
Michael Borisov <br />
Group Id: G

Link to the repository: https://github.com/SecTedd/PSEMolDyn_GroupG <br />
Branch: master <br />
Commit Id: # <br />

Video: <br />

1. Video <br />
   membrane.avi is the simulation of the membrane (task 1) <br />

2. Video <br />
   RT3d.avi shows the Rayleigh-Taylor instability in 3D (task 3) <br />

3. Video <br />
   nano_scale_flow.avi shows the default nano-scale flow simulation until t = 50 (task 4) <br />

## Input
   All input files from assignment 5 can be found in the `/input/` folder <br />

## Build and run

### Build without Doxygen

1. Navigate into the build folder: `cd build`
2. Run `cmake ..` to create the makefile. We use gcc 11.3.0 as compiler.
3. Run `make` to generate the executable and tests.

### Build with Doxygen

1. Navigate into the build folder: `cd build`
2. Run `cmake -D BUILD_DOC=ON ..` to create the makefile. We use gcc 11.3.0 as the compiler.
3. Run `make` to generate the executable and tests and `make doc_doxygen` to generate the documentation.

### Run

After building the project you can run the executable

1. Familiarize yourself with the provided .xml file located in the `/input/` folder and the .xsd file located in the `/xsd/` folder. <br />
   XML-Files are used to specify different parameters. <br />
2. Create a new XML-File or use one of the provided ones and specify all the input parameters that are needed. <br />
3. Navigate into the build folder: `cd build`
4. The application can be run from the console. Run `./MolSim -h` to print the help text. <br />
   The -f option automatically distinguishes between different input files, such as xml, cuboid, sphere etc. <br />
5. The application also allows you to enter an interactive menu where you can read in multiple files, start the simulation multiple times, etc. Run `./MolSim -m` to enter the menu, a help message is shown. <br />
   The -f option automatically distinguishes between different input files, such as xml, cuboid, sphere etc. <br />
6. The generated VTK-Files can be found at the path you specified in the xml. <br />
7. **Warning** the contents of the output folders will be overwritten if you run the simulation multiple times with the same output folder path! <br />

### Checkpoints 

When running the program with the flag `-c` the program creates a checkpoint in the `checkpoints` folder, which has the `PSEMOLDYN_GROUPG` folder as its parent directory. 
If the checkpoint folder does not exist, you need to create it manually as a direct child of the parent directory. 
After writing a checkpoint file, you can read it with the `-f` command or include the path in the filename bracket in a XML file for the configuration. The exact syntax can be found at `src/xsd/Simulation.xsd`

### CSV-Output
When specifying a csvWriteFrequency greater than 0, the program writes .csv files in the `/csv/`directory. <br />
This directory must exist on the same level as the `/src/`folder and needs to be created if it's not there. Otherwise no csv output will be written. <br />

### Run tests

1. Navigate into the build folder: `cd build`
2. Run `ctest` to execute all unit tests.
3. The result of the unit tests are printed to the console.

## Logging

For logging we use spdlog. <br />
The logs are written to files which can be found in the **/logs/** folder or to the console (see `./MolSim -h`). The logs are separated into logic and memory logs. <br />
Logic logs are used to log events in the program flow. Within the logic logs, there is the distinction between input, output and simulation. <br />
Memory logs on the other hand document the construction and destruction of objects and therefore help to detect and prevent memory leaks. <br />

## Contest 2 (Parallel)

Ran with OMP_NUM_THREADS=28 on cm2_tiny and 1 benchmark iteration <br />
GCC: <br />
- 2D 38 sec; 2,61 * 10^6 MUPS <br />
- 3D 1381 sec; 723787 MUPS <br />
ICPX: <br />
- 2D: 42 sec; 2,38 * 10^6 MUPS <br />
- 3D: 1534 sec; 651890 MUPS <br />

## Structure: 
```
./
├── CMakeLists.txt
├── Doxyfile
├── README.md
├── Thread_Benchmark.txt
├── build
├── checkpoints
├── cmake
│   └── modules
│       ├── doxygen.cmake
│       ├── googletest.cmake
│       └── spdlog.cmake
├── csv
├── input
│   ├── Benchmark.xml
│   ├── Membrane.xml
│   ├── Nano-scale-flow.xml
│   ├── RayleighTaylorBig2D.xml
│   ├── RayleighTaylorBig3D.xml
│   ├── Simulation.xml
│   ├── SimulationState.xml
│   ├── eingabe-cuboid.txt
│   ├── eingabe-sonne.txt
│   └── eingabe-sphere.txt
├── libs
├── plots
├── src
│   ├── ConsoleMenu.cpp
│   ├── ConsoleMenu.h
│   ├── MolSim.cpp
│   ├── inputReader
│   │   ├── CuboidInputReader.cpp
│   │   ├── CuboidInputReader.h
│   │   ├── FileReader.cpp
│   │   ├── FileReader.h
│   │   ├── InputFacade.cpp
│   │   ├── InputFacade.h
│   │   ├── InputReader.cpp
│   │   ├── InputReader.h
│   │   ├── SphereInputReader.cpp
│   │   ├── SphereInputReader.h
│   │   ├── XMLInputReader.cpp
│   │   └── XMLInputReader.h
│   ├── model
│   │   ├── Cuboid.cpp
│   │   ├── Cuboid.h
│   │   ├── DirectSumParticleContainer.cpp
│   │   ├── DirectSumParticleContainer.h
│   │   ├── LinkedCellParticleContainer.cpp
│   │   ├── LinkedCellParticleContainer.h
│   │   ├── Particle.cpp
│   │   ├── Particle.h
│   │   ├── ParticleCell.cpp
│   │   ├── ParticleCell.h
│   │   ├── ParticleContainer.h
│   │   ├── ProgramParameters.cpp
│   │   ├── ProgramParameters.h
│   │   ├── Sphere.cpp
│   │   └── Sphere.h
│   ├── outputWriter
│   │   ├── CSVWriter.cpp
│   │   ├── CSVWriter.h
│   │   ├── CheckpointWriter.cpp
│   │   ├── CheckpointWriter.h
│   │   ├── OutputFacade.cpp
│   │   ├── OutputFacade.h
│   │   ├── VTKWriter.cpp
│   │   ├── VTKWriter.h
│   │   ├── XYZWriter.cpp
│   │   ├── XYZWriter.h
│   │   ├── vtk-unstructured.cpp
│   │   ├── vtk-unstructured.h
│   │   └── vtk-unstructured.xsd
│   ├── simulation
│   │   ├── DVProfileCalculator.cpp
│   │   ├── DVProfileCalculator.h
│   │   ├── InterParticleForce.cpp
│   │   ├── InterParticleForce.h
│   │   ├── InterParticleGravitationalForce.cpp
│   │   ├── InterParticleGravitationalForce.h
│   │   ├── LennardJonesForce.cpp
│   │   ├── LennardJonesForce.h
│   │   ├── LennardJonesForceHarmonic.cpp
│   │   ├── LennardJonesForceHarmonic.h
│   │   ├── Simulation.cpp
│   │   ├── Simulation.h
│   │   ├── SingleParticleForce.cpp
│   │   ├── SingleParticleForce.h
│   │   ├── SingleParticleGravitationalForce.cpp
│   │   ├── SingleParticleGravitationalForce.h
│   │   ├── TemporalSingleParticleForce.cpp
│   │   ├── TemporalSingleParticleForce.h
│   │   ├── Thermostat.cpp
│   │   └── Thermostat.h
│   ├── utils
│   │   ├── ArrayUtils.h
│   │   ├── File.h
│   │   ├── Input.h
│   │   ├── Logger.h
│   │   ├── MaxwellBoltzmannDistribution.h
│   │   ├── PContainer.h
│   │   └── ParticleGenerator.h
│   └── xsd
│       ├── Simulation.cxx
│       ├── Simulation.hxx
│       ├── Simulation.xsd
│       ├── SimulationState.cxx
│       ├── SimulationState.hxx
│       └── SimulationState.xsd
└── tests
    ├── CSVWriter_test.cc
    ├── Checkpoints_test.cc
    ├── CuboidInputReader_test.cc
    ├── DVCalculator_test.cc
    ├── DirectSumParticleContainer_test.cc
    ├── GetNeighbours_test.cc
    ├── LennardJonesForce_test.cc
    ├── LinkedCellParticleContainer_test.cc
    ├── Membrane_test.cc
    ├── Parallelization_test.cc
    ├── ParticleCell_test.cc
    ├── SigmaEpsilon_test.cc
    ├── Simulation.xml
    ├── SimulationState.xml
    ├── SingleParticleGravitationalForce_test.cc
    ├── SphereInputReader_test.cc
    ├── StationaryParticles_test.cc
    ├── Thermostat_test.cc
    ├── XMLInputReader_test.cc
    ├── checkpoints
    ├── csv
    ├── eingabe-cuboid.txt
    ├── eingabe-sphere.txt
    └── main.cc
```
