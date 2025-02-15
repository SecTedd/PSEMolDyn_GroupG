#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "../src/model/DirectSumParticleContainer.h"
#include "../src/inputReader/XMLInputReader.h"

// correctness of xml reading for simulation object
TEST(XMLInputReader, XMLSimulation)
{
    std::shared_ptr<ProgramParameters> pp = std::make_shared<ProgramParameters>();
    std::unique_ptr<XMLInputReader> xml = std::make_unique<XMLInputReader>();

    const char *file = "../tests/Simulation.xml";
    xml->readInput(*pp, file);

    EXPECT_THAT(pp->getEndTime(), 10);
    EXPECT_THAT(pp->getDeltaT(), 0.00005);
    EXPECT_THAT(pp->getCutoff(), 3.0);
    EXPECT_THAT(pp->getDimension(), 3);
    BoundaryCondition b = pp->getBoundaries()[2];
    EXPECT_THAT(b, BoundaryCondition::Reflecting);
    EXPECT_THAT(pp->getTempInit(), 0);
    EXPECT_THAT(pp->getBrownianMotion(), false);
    EXPECT_THAT(pp->getNThermostats(), 0);
    EXPECT_THAT(pp->getTempTarget(), 0);
    EXPECT_THAT(pp->getMembrane(), true);
    EXPECT_THAT(pp->getCsvWriteFrequency(), 10000);
    EXPECT_THAT(pp->getNumBins(), 30);
    EXPECT_THAT(pp->getForces().size(), 1);
}