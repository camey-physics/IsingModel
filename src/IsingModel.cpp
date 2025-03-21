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

double IsingModel::getBeta() const{
    return beta_;
}

void IsingModel::setSpin(int i, int j, int k, int val) {
    if (val != 1 && val != -1) {
        throw std::invalid_argument("Spin value must be +1 or -1");
    }
    int local_h = calcLocalH(index(i, j, k));
    // double deltaE = calcDeltaE(index(i, j, k));
    spins_[index(i, j, k)] = val;
    energy_ += -2 *J_ *local_h *val; // deltaE = -2*J*s_i *sum_j(s_<ij>)
}

void IsingModel::setBeta(double beta) {
    beta_ = beta;
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
        for (int n = 0; n < 6; n += 2) {
            int j = NT_[i *6 + n];
            energy_ += spins_[i] *spins_[j];
        }
    }
    energy_ *= -J_;
}

int IsingModel::calcLocalH(int i) const {
    int local_h = 0;
    for (int n = 0; n < 6; ++n) {
        int j = NT_[i *6 + n];
        local_h += spins_[j];
    }
    return local_h;
}

void IsingModel::metropolis(int i) {
    // Can pre-calculate and save a weight table so we don't recalculate each time
    double deltaE = 2 * J_ * spins_[i] * calcLocalH(i);
    if (deltaE <= 0 || gsl_rng_uniform(r) < exp(-beta_ * deltaE)) {
        spins_[i] *= -1;
    }
}

void IsingModel::heatBath(int i) {
    int localH = calcLocalH(i);
    double probUp = 1 /(1 + exp(-2 *beta_ *J_ *localH));
    if (gsl_rng_uniform(r) < probUp) {
        spins_[i] = 1;
    }
    else {
        spins_[i] = -1;
    }
}

void IsingModel::MonteCarloSweep(int numSweeps, bool sequential, void (IsingModel::*update)(int)) {
    if (sequential) {
        // Sequential update
        for (int sweep = 0; sweep < numSweeps; ++sweep) {
            for (int s = 0; s < L_ * L_ * L_; ++s) {
                (this->*update)(s);
            }
        }
    } else {
        // Random update
        for (int sweep = 0; sweep < numSweeps; ++sweep) {
            for (int s = 0; s < L_ * L_ * L_; ++s) {
                int ind = gsl_rng_uniform_int(r, L_ * L_ * L_);
                (this->*update)(ind);
            }
        }
    }
}
