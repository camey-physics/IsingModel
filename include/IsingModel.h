#ifndef ISING_MODEL_H
#define ISING_MODEL_H

#include <vector>
#include <gsl/gsl_rng.h>
#include <memory>

class IsingModel {
public:
    explicit IsingModel(int L = 4, double beta = 1, int seed = 5000, double J = 1); // Constructor
    ~IsingModel();

    int getSpin(int i, int j, int k) const; // Accessor for spins
    std::vector<int> getNeighbors(int index) const;
    double getEnergy() const;
    double getBeta() const;

    void setSpin(int i, int j, int k, int val);
    void setBeta(double beta);
    void calcEnergy();
    void metropolis(int i);
    void heatBath(int i);
    void MonteCarloSweep(int numSweeps, bool sequential=0, void (IsingModel::*update)(int) = &IsingModel::metropolis);

private:
    int L_;
    double J_;
    double beta_;
    int seed_;
    std::vector<int> spins_;
    std::vector<int> NT_;

    double energy_;

    inline int index(int i, int j, int k) const;
    inline int mod(int i) const;
    void initializeNT();
    gsl_rng *r;
    int calcLocalH(int i) const;

};

#endif // ISING_MODEL_H
