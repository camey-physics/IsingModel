#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_rng.h>
#include "../include/IsingModel.h"

#include <vector>
#include <gsl/gsl_rng.h>
#include <iostream>
#include <cmath>

std::pair<double, double> bootstrapBinderCumulant(
    const std::vector<double>& m2, 
    const std::vector<double>& m4, 
    gsl_rng* r, 
    int numBootstraps = 1000) 
{
    int N = m2.size();
    if (N != m4.size() || N == 0) {
        throw std::invalid_argument("Vectors m2 and m4 must have the same nonzero size.");
    }

    std::vector<double> bootCumulants(numBootstraps);

    // Perform bootstrap resampling
    for (int b = 0; b < numBootstraps; ++b) {
        double sumM2 = 0.0, sumM4 = 0.0;
        for (int i = 0; i < N; ++i) {
            int idx = gsl_rng_uniform_int(r, N); // Random index in range [0, N-1]
            sumM2 += m2[idx];
            sumM4 += m4[idx];
        }
        double avgM2 = sumM2 / N;
        double avgM4 = sumM4 / N;
        bootCumulants[b] = 1.0 - (avgM4 / (3.0 * avgM2 * avgM2));
    }

    // Compute mean of the bootstrap cumulants
    double meanCumulant = 0.0;
    for (int b = 0; b < numBootstraps; ++b) {
        meanCumulant += bootCumulants[b];
    }
    meanCumulant /= numBootstraps;

    // Compute standard deviation (bootstrap error)
    double variance = 0.0;
    for (int b = 0; b < numBootstraps; ++b) {
        double diff = bootCumulants[b] - meanCumulant;
        variance += diff * diff;
    }
    double stdError = std::sqrt(variance / (numBootstraps - 1));

    return {meanCumulant, stdError};
}


void runSimulatedAnnealing(IsingModel &model, double betaStart, double betaEnd, double betaStep, int numSweeps, int equilibrationSweeps, int numMeasurements, int seed = 6000) {
    std::ofstream outFile("binder_cumulant_data.txt");
    outFile << "# Beta U_L\n";  // Column headers for beta and Binder cumulant
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, seed);

    std::vector<double>  m2, m4;
    m2.resize(numMeasurements,0);
    m4.resize(numMeasurements,0);
    
    // Run annealing from betaStart to betaEnd
    for (double beta = betaStart; beta <= betaEnd; beta += betaStep) {
        // Equilibrate the system at current beta
        model.setBeta(beta);
        model.monteCarloSweep(equilibrationSweeps, true, &IsingModel::heatBath); // Equilibrate

        // Take measurements at current beta
        for (int i = 0; i < numMeasurements; ++i) {
            // Perform a Monte Carlo sweep to update the system
            model.monteCarloSweep(numSweeps, true, &IsingModel::heatBath); // Sweep

            double magnetization = model.calcMagnetization();
            m2[i] = magnetization *magnetization;
            m4[i] = m2[i] *m2[i];
        }

        // Calculate the Binder cumulant
        auto [binderCumulant, binderError] = bootstrapBinderCumulant(m2, m4, r, 6000);

        // Output the results for this beta
        outFile << beta << " " << binderCumulant << " " << binderError << std::endl;
        std::cout << "Beta: " << beta << ", Binder Cumulant: " << binderCumulant 
                  << ", Error: " << binderError << std::endl;
        
    }

    outFile.close();
}

int main() {
    // Initialize Ising model
    IsingModel model(8, 1.0, 5000, 1.0);

    // Simulated annealing parameters
    double betaStart = 0.15;
    double betaEnd = 0.25;
    double betaStep = 0.01;
    int numSweeps = 1000;
    int equilibrationSweeps = 100000;
    int numMeasurements = 1000;

    // Run the simulated annealing and measurement of Binder cumulant
    runSimulatedAnnealing(model, betaStart, betaEnd, betaStep, numSweeps, equilibrationSweeps, numMeasurements);

    return 0;
}
