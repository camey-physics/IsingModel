#ifndef ISING_MODEL_H
#define ISING_MODEL_H

#include <vector>
#include <gsl/gsl_rng.h>

struct IsingParams {
    int L;
    double J;
    double beta;
    int seed;
    
    bool operator==(const IsingParams& other) const {
        return L == other.L && J == other.J && seed == other.seed;
    }
};

class IsingModel {
public:
    explicit IsingModel(int L = 4, double J = 1, double beta = 1, int seed = 5000); // Constructor
    ~IsingModel();

    IsingParams getParameters() const;
    int getSpin(int i, int j, int k) const; // Accessor for spins
    std::vector<int> getNeighbors(int index) const;
    double getEnergy() const;

    void setSpin(int i, int j, int k, int val);

private:
    IsingParams params_;
    int L_;
    double J_;
    double beta_;
    int seed_;
    std::vector<int> spins_;
    std::vector<int> NT_;

    double energy_;

    int index(int i, int j, int k) const;
    int mod(int i) const;
    void initializeNT();
    gsl_rng *r;
    void calcEnergy();
    double calcDeltaE(int i);
    void metropolis(int i);

};

#endif // ISING_MODEL_H
