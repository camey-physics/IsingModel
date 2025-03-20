#include "../include/IsingModel.h"
#include <gsl/gsl_rng.h>
#include <stdexcept>
#include <cmath>

IsingModel::IsingModel(int L, double J, double beta, int seed) : L_(L), J_(J), beta_(beta), seed_(seed), spins_(L *L *L, 1), params_{L, J, beta, seed} {
    r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, seed_);
    // for (int i = 0; i < L_ * L_ * L_; ++i) {
    //     spins_[i] = (gsl_rng_uniform(r) < 0.5) ? -1 : 1;
    // }
    // energy_ = -J_ *L_ *L_ *L_;
    initializeNT();
    calcEnergy();
}

IsingModel::~IsingModel() {
    gsl_rng_free(r);
}

IsingParams IsingModel::getParameters() const {
    return params_;
}

int IsingModel::getSpin(int i, int j, int k) const {
    return spins_[index(i, j, k)];
}

int IsingModel::index(int i, int j, int k) const {
    return mod(i) * L_ * L_ + mod(j) * L_ + mod(k);
}
std::vector<int> IsingModel::getNeighbors(int ind) const { // Only used for testing NT_. For actual runs, directly calculate indices for speed
    std::vector<int> neighbors(6);
    for (int i = 0; i < 6; ++i) {
        neighbors[i] = NT_[ind*6+i];
    }
    return neighbors;
}

double IsingModel::getEnergy() const {
    return energy_;
}

void IsingModel::setSpin(int i, int j, int k, int val) {
    if (val != 1 && val != -1) {
        throw std::invalid_argument("Spin value must be +1 or -1");
    }
    double deltaE = calcDeltaE(index(i, j, k));
    spins_[index(i, j, k)] = val;
    energy_ += deltaE;
}

int IsingModel::mod(int i) const {
        return (i % L_ + L_) % L_;   
}

void IsingModel::initializeNT() {
    NT_.resize(L_ *L_ *L_ *6);
    for (int i = 0; i < L_; ++i) {
        for (int j = 0; j < L_; ++j) {
            for (int k = 0; k < L_; ++k) {
                int ind = index(i, j, k);
                NT_[ind*6+0] = index(i-1, j, k);
                NT_[ind*6+1] = index(i+1, j, k);
                NT_[ind*6+2] = index(i, j-1, k);
                NT_[ind*6+3] = index(i, j+1, k);
                NT_[ind*6+4] = index(i, j, k-1);
                NT_[ind*6+5] = index(i, j, k+1);
            }
        }
    }
}

void IsingModel::calcEnergy() {
    energy_ = 0.0;
    for (int i = 0; i < L_ *L_ *L_; ++i) {
        for (int n = 0; n < 6; ++n) {
            int j = NT_[i *6 + n];
            energy_ += spins_[i] *spins_[j];
        }
    }
    energy_ *= -J_ /2.0;
}


double IsingModel::calcDeltaE(int i) {
    double deltaE = 0.0;
    for (int n = 0; n < 6; ++n) {
        int j = NT_[i *6 + n];
        deltaE += spins_[j];
    }
    deltaE *= 2 *J_ *spins_[i];
    return deltaE;
}

void IsingModel::metropolis(int i) {
    double deltaE = calcDeltaE(i); // Can calculate and save a weight table so we don't recalculate weights
    double weight = exp(-beta_ *deltaE);
    if (weight >= 1.0 || gsl_rng_uniform(r) < weight) {
        spins_[i] *= -1;
    }
}