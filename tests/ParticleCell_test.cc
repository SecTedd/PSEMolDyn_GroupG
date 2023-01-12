#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "../src/model/ParticleCell.h"
#include "../src/utils/ArrayUtils.h"

/**
 * test replacement of invalid particles when new particle is inserted
 */
TEST(ParticleCell, InsertParticle)
{
    BoundaryCondition outflow = BoundaryCondition::Outflow;
    std::array<BoundaryCondition, 6> boundaries = {outflow, outflow, outflow, outflow, outflow, outflow};
    ParticleCell cell = ParticleCell(CellType::InnerCell, boundaries);

    cell.reserveMemory(3);
    Particle p1 = Particle();
    Particle p2 = Particle();
    Particle p3 = Particle();
    std::vector<Particle> particles;
    particles.push_back(p1);
    particles.push_back(p2);
    particles.push_back(p3);

    cell.insertParticleIndex(0);
    cell.insertParticleIndex(1);

    EXPECT_EQ(cell.size(), 2); // particles should be added normally

    particles[1].setInvalid(true);
    cell.removeInvalid(&particles);
    EXPECT_EQ(cell.size(), 1); // invalid particles are erased

    cell.insertParticleIndex(2);

    EXPECT_EQ(cell.size(), 2); // invalid particle p2 should be overwritten with p3
}