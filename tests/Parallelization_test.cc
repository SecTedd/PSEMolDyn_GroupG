#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "../src/model/LinkedCellParticleContainer.h"

TEST(CellGrouping, NoParallelization) {
    double cutoff = 3;
    std::array<double, 3> domain = {6, 6, 6};
    BoundaryCondition out = BoundaryCondition::Outflow;
    std::array<BoundaryCondition, 6> boundaries = {out, out, out, out, out, out};
    int parallelization = 0; //since tests are compiled without OpenMP this would automatically reset to 0 anyways
    LinkedCellParticleContainer pc = LinkedCellParticleContainer(cutoff, domain, boundaries, parallelization);
    std::vector<std::vector<int>> &cellGroups = pc.getCellGroups();
    EXPECT_EQ(pc.getParallel(), 0);
    EXPECT_EQ(cellGroups.size(), 1);
    EXPECT_EQ(cellGroups[0].size(), 8);
    EXPECT_THAT(cellGroups[0], testing::ElementsAre(21, 22, 25, 26, 37, 38, 41, 42));
}

TEST(CellGrouping, ParallelizationStrategy1) {
    double cutoff = 3;
    std::array<double, 3> domain = {9, 9, 9};
    BoundaryCondition out = BoundaryCondition::Outflow;
    std::array<BoundaryCondition, 6> boundaries = {out, out, out, out, out, out};
    int parallelization = 0; //since tests are compiled without OpenMP this would automatically reset to 0 anyways
    LinkedCellParticleContainer pc = LinkedCellParticleContainer(cutoff, domain, boundaries, parallelization);
    //reset parallelization strategy to 1 and reinitialize cells for correct grouping
    pc.setParallel(1);
    pc.initializeCells(boundaries, 1);
    std::vector<std::vector<int>> &cellGroups = pc.getCellGroups();
    EXPECT_EQ(pc.getParallel(), 1);
    EXPECT_EQ(cellGroups.size(), 18);
    for (auto &group : cellGroups) {
        EXPECT_EQ(group.capacity(), 2);
    }
    EXPECT_THAT(cellGroups[0], testing::ElementsAre(31, 81));
    EXPECT_THAT(cellGroups[1], testing::ElementsAre(32, 82));
    EXPECT_THAT(cellGroups[2], testing::ElementsAre(33, 83));
    EXPECT_THAT(cellGroups[3], testing::ElementsAre(36, 86));
    EXPECT_THAT(cellGroups[4], testing::ElementsAre(37, 87));
    EXPECT_THAT(cellGroups[5], testing::ElementsAre(38, 88));
    EXPECT_THAT(cellGroups[6], testing::ElementsAre(41, 91));
    EXPECT_THAT(cellGroups[7], testing::ElementsAre(42, 92));
    EXPECT_THAT(cellGroups[8], testing::ElementsAre(43, 93));
    EXPECT_THAT(cellGroups[9], testing::ElementsAre(56));
    EXPECT_THAT(cellGroups[10], testing::ElementsAre(57));
    EXPECT_THAT(cellGroups[11], testing::ElementsAre(58));
    EXPECT_THAT(cellGroups[12], testing::ElementsAre(61));
    EXPECT_THAT(cellGroups[13], testing::ElementsAre(62));
    EXPECT_THAT(cellGroups[14], testing::ElementsAre(63));
    EXPECT_THAT(cellGroups[15], testing::ElementsAre(66));
    EXPECT_THAT(cellGroups[16], testing::ElementsAre(67));
    EXPECT_THAT(cellGroups[17], testing::ElementsAre(68));
}
