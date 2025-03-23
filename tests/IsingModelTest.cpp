#include "../include/IsingModel.h"
#include <gtest/gtest.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <tuple>

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

TEST(IsingModelTest, TestParams) {
    IsingModel model;
    auto [L, J, beta, seed] = model.getParams();
    EXPECT_EQ(L, 4);
    EXPECT_EQ(beta, 1.0);
    EXPECT_EQ(J, 1.0);
    EXPECT_EQ(seed, 5000);
}

TEST(IsingModelTest, CalcEnergy) {
    int L = 5;
    IsingModel model(L);
    EXPECT_NEAR(model.calcEnergy(), -3.0, 1e-10);
    model.setSpin(1,0,0,-1);
    EXPECT_NEAR(model.calcEnergy(), -3.0 + 12.0 /L /L /L, 1e-10);
    model.setSpin(0,0,0,-1);
    EXPECT_NEAR(model.calcEnergy(), -3.0 + 20.0 /L /L /L, 1e-10);
    model.setSpin(0,L-1,0,-1);
    EXPECT_NEAR(model.calcEnergy(), -3.0 + 28.0 /L /L /L, 1e-10);
}

TEST(IsingModelTest, MetropolisSweep) {
    int L = 10;
    int numSamples = 100;
    double beta = 0.1;
    IsingModel model(L);
    model.setBeta(beta);
    double avg_energy = 0.0, avg_mag = 0.0;
    for (int i = 0; i < numSamples; ++i) {
        model.monteCarloSweep(100);
        avg_energy += model.calcEnergy();
        avg_mag += model.calcMagnetization();
    }
    avg_energy /= numSamples;
    avg_mag /= numSamples;
    
    EXPECT_NEAR(avg_energy, -3/2 *tanh(3*beta), 5e-2);
    EXPECT_NEAR(avg_mag, 0.0, 5e-2);
}

TEST(IsingModelTest, HeatBathSweep) {
    int L = 10;
    int numSamples = 100;
    double beta = 0.1;
    IsingModel model(L);
    model.setBeta(beta);
    double avg_energy = 0.0, avg_mag = 0.0;
    for (int i = 0; i < numSamples; ++i) {
        model.monteCarloSweep(100, 0, &IsingModel::heatBath);
        avg_energy += model.calcEnergy();
        avg_mag += model.calcMagnetization();
    }
    avg_energy /= numSamples;
    avg_mag /= numSamples;
    
    EXPECT_NEAR(avg_energy, -3/2 *tanh(3*beta), 5e-2);
    EXPECT_NEAR(avg_mag, 0.0, 5e-2);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
