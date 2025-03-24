#include <iostream>
#include <fstream>
#include <cmath>
#include "../include/IsingModel.h"

void runSimulatedAnnealing(IsingModel &model, double betaStart, double betaEnd, double betaStep, int numSweeps, int equilibrationSweeps, int numMeasurements) {
    std::ofstream outFile("binder_cumulant_data.txt");
    outFile << "# Beta U_L\n";  // Column headers for beta and Binder cumulant

    // Run annealing from betaStart to betaEnd
    for (double beta = betaStart; beta <= betaEnd; beta += betaStep) {
        // Equilibrate the system at current beta
        model.setBeta(beta);
        model.monteCarloSweep(equilibrationSweeps, true, &IsingModel::heatBath); // Equilibrate

        double sumM2 = 0.0, sumM4 = 0.0;

        // Take measurements at current beta
        for (int i = 0; i < numMeasurements; ++i) {
            // Perform a Monte Carlo sweep to update the system
            model.monteCarloSweep(numSweeps, true, &IsingModel::heatBath); // Sweep

            double magnetization = model.calcMagnetization();
            double m2 = magnetization * magnetization;
            double m4 = m2 * m2;

            sumM2 += m2;
            sumM4 += m4;
        }

        // Calculate the Binder cumulant
        double binderCumulant = 1.0 - (sumM4 / numMeasurements) / (3.0 * std::pow(sumM2 / numMeasurements, 2));

        // Output the results for this beta
        outFile << beta << " " << binderCumulant << std::endl;
        std::cout << "Beta: " << beta << ", Binder Cumulant: " << binderCumulant << std::endl;
    }

    outFile.close();
}

int main() {
    // Initialize Ising model
    IsingModel model(5, 1.0, 5000, 1.0);

    // Simulated annealing parameters
    double betaStart = 0.05;
    double betaEnd = 0.3;
    double betaStep = 0.025;
    int numSweeps = 10000;
    int equilibrationSweeps = 5000;
    int numMeasurements = 100;

    // Run the simulated annealing and measurement of Binder cumulant
    runSimulatedAnnealing(model, betaStart, betaEnd, betaStep, numSweeps, equilibrationSweeps, numMeasurements);

    return 0;
}
