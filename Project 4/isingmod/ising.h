#ifndef ISING_H
#define ISING_H

//Function for Metropolis algo, also loops on MC cycles
void MetropolisSampling(int, int, double, double*);

//Define the lattice
void InitializeLattice(int, double**, double&, double&);

//Print to file function
void WriteResultstoFile(int, int, double, double*,int);

//inline for PBCs
inline int PeriodicBoundary(int i, int limit, int add) {
  return (i+limit+add) % (limit);
}
void Update(double* ExpectationValues, double Energy, double MagneticMoment);

void PrintEnergyCounts(int EnergyCounts[], double Temperature);

#endif // ISING_H
