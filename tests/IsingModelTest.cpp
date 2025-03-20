#include "../include/IsingModel.h"
#include <gtest/gtest.h>
#include <vector>

TEST(IsingModelTest, InitializesSpinLattice) {
    int L = 4;
    IsingModel model(L);

    for (int i = 0; i < L-1; ++i) {
        for (int j = 0; j < L-1; ++j) {
            for (int k = 0; k < L-1; ++k) {
                EXPECT_TRUE(model.getSpin(i, j, k) == +1 || model.getSpin(i,j,k) == -1);
            }
        }
    }
}

TEST(IsingModelTest, PeriodicBoundaryConditions) {
    int L = 4;
    IsingModel model(L);
    // for (int i = 0; i < L; ++i) {
    //     for (int j = 0; j < L; ++j) {
    //         for (int k = 0; k < L; ++k) {
    //             model.setSpin(i, j, k, 1);
    //         }
    //     }
    // }
    model.setSpin(L-1, L-1, L-1, -1);
    EXPECT_TRUE(model.getSpin(-1, -1, -1) == model.getSpin(L-1, L-1, L-1));
}

TEST(IsingModelTest, NeighborTable) {
    int L = 4;
    IsingModel model(L);
    std::vector<int> expectedNeighbors = {(L-1)*L*L, L*L, (L-1)*L, L, L-1, 1};
    EXPECT_EQ(model.getNeighbors(0), expectedNeighbors);
    expectedNeighbors = {0, 32, 28, 20, 19, 17};
    EXPECT_EQ(model.getNeighbors(L*L), expectedNeighbors);
}

TEST(IsingModelTest, DefaultParameters) {
    IsingModel model;
    IsingParams defaultParams{4, 1.0, 1.0, 5000};
    EXPECT_EQ(model.getParameters(), defaultParams);
}

TEST(IsingModelTest, CalcEnergy) {
    int L = 4;
    IsingModel model(L);
    EXPECT_EQ(model.getEnergy(), -3 *L *L *L);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
