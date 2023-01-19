#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "../src/model/LinkedCellParticleContainer.h"
#include "../src/inputReader/XMLInputReader.h"
#include "../src/outputWriter/CheckpointWriter.h"
#include "../src/utils/File.h"

// correctness of xml reading for simulation state
TEST(Checkpoints, ReadCheckpoint)
{
    std::shared_ptr<ProgramParameters> pp = std::make_shared<ProgramParameters>();
    std::unique_ptr<XMLInputReader> xml = std::make_unique<XMLInputReader>();

    const char *file = "../tests/SimulationState.xml";
    xml->readInput(*pp, file);
    Particle p = pp->getParticleContainer()->getActiveParticles()[0];
    EXPECT_THAT(p.getM(), 1);
    EXPECT_THAT(p.getEpsilon(), 2);
    EXPECT_THAT(p.getSigma(), 3);
    EXPECT_THAT(p.getType(), 4);
    std::array<double, 3> x;
    x = p.getX();
    EXPECT_THAT(x[0], 1);
    std::array<double, 3> v;
    v = p.getV();
    EXPECT_THAT(v[0], 2);
    std::array<double, 3> f;
    f = p.getF();
    EXPECT_THAT(f[0], 3);
    std::array<double, 3> old_f;
    old_f = p.getOldF();
    EXPECT_THAT(old_f[0], 4);
}

// correctness of xml writing and reading for simulation state
TEST(Checkpoints, WriteAndReadCheckpoint)
{
    std::shared_ptr<ProgramParameters> pp = std::make_shared<ProgramParameters>();
    std::unique_ptr<XMLInputReader> xml = std::make_unique<XMLInputReader>();
    std::unique_ptr<outputWriter::CheckpointWriter> cp = std::make_unique<outputWriter::CheckpointWriter>();
    std::string filename = "../tests/checkpoints/" + File::getDateTime() + "test_checkpoint.xml";
    std::string scheme = "../../src/xsd/SimulationState.xsd"; 

    double cutoff = 3;
    std::array<double, 3> domain = {9, 9, 9};
    BoundaryCondition out = BoundaryCondition::Outflow;
    std::array<BoundaryCondition, 6> boundaries = {out, out, out, out, out, out};
    std::unique_ptr<ParticleContainer> pc;
    pc.reset(new LinkedCellParticleContainer(cutoff, domain, boundaries));

    std::array<double, 3> x = {1, 2, 0};
    std::array<double, 3> v = {-1, -2, -3};
    std::array<double, 3> f = {-4, -5, -6};
    std::array<double, 3> oldF = {-7, -8, -9};
    double mass = 4;
    double sigma = 5;
    double epsilon = 6;
    int type = 7;
    double stiffness = 8; 
    double averageBondLength = 9; 

    pc->addParticle(x, v, f, oldF, mass, epsilon, sigma, type, stiffness, averageBondLength);
    cp->writeCheckpoint(pc.get(), &filename, &scheme);

    xml->readInput(*pp, filename.c_str());
    Particle p = pp->getParticleContainer()->getActiveParticles()[0];
    EXPECT_THAT(p.getM(), 4);
    EXPECT_THAT(p.getEpsilon(), 6);
    EXPECT_THAT(p.getSigma(), 5);
    EXPECT_THAT(p.getType(), 7);
    EXPECT_THAT(p.getStiffness(), 8); 
    EXPECT_THAT(p.getAverageBondLength(), 9); 
    EXPECT_THAT(p.getX(), testing::ElementsAre(1, 2, 0));
    EXPECT_THAT(p.getV(), testing::ElementsAre(-1, -2, -3));
    EXPECT_THAT(p.getF(), testing::ElementsAre(-4, -5, -6));
    EXPECT_THAT(p.getOldF(), testing::ElementsAre(-7, -8, -9));
}