#include "../include/IsingModel.h"
#include <gtest/gtest.h>
#include <vector>
#include <iostream>
#include <cmath>

TEST(IsingModelTest, InitializesSpinLattice) {
    int L = 4;
    IsingModel model(L);

    for (int i = 0; i < L; ++i) {
        for (int j = 0; j < L; ++j) {
            for (int k = 0; k < L; ++k) {
                EXPECT_TRUE(model.getSpin(i, j, k) == +1 || model.getSpin(i,j,k) == -1);
            }
        }
    }
}

TEST(IsingModelTest, PeriodicBoundaryConditions) {
    int L = 4;
    IsingModel model(L);
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

TEST(IsingModelTest, CalcEnergy) {
    int L = 5;
    IsingModel model(L);
    EXPECT_EQ(model.getEnergy(), -3 *L *L *L);
    model.setSpin(1,0,0,-1);
    EXPECT_EQ(model.getEnergy(), -3 *L *L *L + 12);
    model.setSpin(0,0,0,-1);
    EXPECT_EQ(model.getEnergy(), -3 *L *L *L + 20);
    model.setSpin(0,L-1,0,-1);
    EXPECT_EQ(model.getEnergy(), -3 *L *L *L + 28);
}

TEST(IsingModelTest, MetropolisSweep) {
    int L = 10;
    int numSamples = 100;
    double beta = 0.1;
    IsingModel model(L);
    model.setBeta(beta);
    double avg_energy = 0.0;
    for (int i = 0; i < numSamples; ++i) {
        model.MonteCarloSweep(100);
        model.calcEnergy();
        avg_energy += model.getEnergy() /L /L /L;
    }
    avg_energy /= numSamples;
    
    EXPECT_NEAR(avg_energy, -3/2 *tanh(3*beta), 5e-2);
}

TEST(IsingModelTest, HeatBathSweep) {
    int L = 10;
    int numSamples = 100;
    double beta = 0.1;
    IsingModel model(L);
    model.setBeta(beta);
    double avg_energy = 0.0;
    for (int i = 0; i < numSamples; ++i) {
        model.MonteCarloSweep(100, 0, &IsingModel::heatBath);
        model.calcEnergy();
        avg_energy += model.getEnergy() /L /L /L;
    }
    avg_energy /= numSamples;
    
    EXPECT_NEAR(avg_energy, -3/2 *tanh(3*beta), 5e-2);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
