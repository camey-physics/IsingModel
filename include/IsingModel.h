#ifndef ISING_MODEL_H
#define ISING_MODEL_H

#include <vector>
#include <gsl/gsl_rng.h>
#include <memory>
#include <tuple>

class IsingModel {
public:
    explicit IsingModel(int L = 4, double beta = 1, int seed = 5000, double J = 1); // Constructor
    ~IsingModel();

    enum class UpdateMethod { metropolis, heatBath, wolff };

    int getSpin(int i, int j, int k) const; // Accessor for spins
    std::vector<int> getNeighbors(int index) const;
    double calcEnergy() const;
    double getBeta() const;
    double calcMagnetization() const;
    auto getParams() const {
        return std::make_tuple(L_, J_, beta_, seed_);
    }
    void setSpin(int i, int j, int k, int val);
    void setBeta(double beta);
    void monteCarloSweep(int numSweeps, UpdateMethod method, bool sequential=0);

private:
    int L_;
    double J_;
    double beta_;
    int seed_;
    std::vector<int> spins_;
    std::vector<int> NT_;

    void metropolis(int i);
    void heatBath(int i);
    int wolff();

    inline int index(int i, int j, int k) const;
    inline int mod(int i) const;
    void initializeNT();
    gsl_rng *r;
    int calcLocalH(int i) const;

};

#endif // ISING_MODEL_H
