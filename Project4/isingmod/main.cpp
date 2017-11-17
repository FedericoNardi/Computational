#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <string>
#include "mpi.h"
using namespace  std;
// output file
ofstream ofile;

// inline function for PeriodicBoundary boundary conditions
inline int PeriodicBoundary(int i, int limit, int add) {
  return (i+limit+add) % (limit);
}
// Function to initialise energy and magnetization
void InitializeLattice(int, double**, double&, double&);
// The metropolis algorithm including the loop over Monte Carlo cycles
void MetropolisSampling(int, int, double, double*);
// prints to file the results of the calculations
void WriteResultstoFile(int, int, double, double*,int);

// Main program begins here

int main(int argc, char* argv[])
{
    string filename = "IsingRandTrange_";
    int NSpins=100;
    int MCcycles=1e6;
    double InitialTemp=2.0;   //kT/J
    double FinalTemp=2.5;
    double TempStep=0.001;
    int numprocs, my_rank;
    // Declare new file name and add lattice size to file name
    //cout<<"MCcycles:" <<MCcycles <<"\n";

    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

    if(my_rank==0){
        string fileout = filename;
        string argument = to_string(NSpins);
        fileout.append(argument);
        fileout.append(".txt");
        ofile.open(fileout);
    }
    // Start Monte Carlo sampling by looping over the selected Temperatures
    for (double Temperature = InitialTemp; Temperature <= FinalTemp; Temperature+=TempStep){
        double* ExpectationValues = new double [6];
        for(int ii=0; ii<=5; ii++) ExpectationValues[ii]=0.;

        // Start Monte Carlo computation and get expectation values
        MetropolisSampling(NSpins, MCcycles, Temperature, ExpectationValues);

        double* TotalExpectationValues = new double [6];
        for(int ii=0; ii<=5; ii++){
            TotalExpectationValues[ii]=0.;
            MPI_Reduce(&ExpectationValues[ii], &TotalExpectationValues[ii], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }

        WriteResultstoFile(NSpins, MCcycles*numprocs, Temperature, TotalExpectationValues, 0);
        delete [] ExpectationValues;
    }
    if(my_rank==0) ofile.close();  // close output file

    MPI_Finalize();
    return 0;
}



// The Monte Carlo part with the Metropolis algo with sweeps over the lattice
void MetropolisSampling(int NSpins, int MCcycles, double Temperature, double* ExpectationValues)
{
    // Initialize the seed and call the Mersenne algo
    std::random_device rd;
    std::mt19937_64 gen(rd());
    // Set up the uniform distribution for x \in [[0, 1]
    std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
    // Initialize the lattice spin values
    double** SpinMatrix = new double*[NSpins];
    for(int i = 0; i<NSpins; i++){
        SpinMatrix[i]=new double[NSpins];
    }
    //    initialize energy and magnetization
    double Energy = 0.;     double MagneticMoment = 0.;
    // initialize array for expectation values
    InitializeLattice(NSpins, SpinMatrix, Energy, MagneticMoment);
    // setup array for possible energy changes
    double* EnergyDifference = new double[5];
    for( int de =-8; de <= 8; de+=4){
        double index = 0.25*(de+8);
        EnergyDifference[(int)index] = exp(-de/Temperature);
    }

    // Start Monte Carlo cycles
    for (int cycles = 1; cycles <= MCcycles; cycles++){
        int AcceptedCycles;
        // The sweep over the lattice, looping over all spin sites
        for(int x =0; x < NSpins; x++) {
            for (int y= 0; y < NSpins; y++){
                int ix = (int) (RandomNumberGenerator(gen)*(double)NSpins);
                int iy = (int) (RandomNumberGenerator(gen)*(double)NSpins);
                int deltaE =  2*SpinMatrix[ix][iy]*
                        (SpinMatrix[ix][PeriodicBoundary(iy,NSpins,-1)]+
                        SpinMatrix[PeriodicBoundary(ix,NSpins,-1)][iy] +
                        SpinMatrix[ix][PeriodicBoundary(iy,NSpins,1)] +
                        SpinMatrix[PeriodicBoundary(ix,NSpins,1)][iy]);
                if (EnergyDifference<=0 || RandomNumberGenerator(gen)<=EnergyDifference[(deltaE+8)/4]) {
                    SpinMatrix[ix][iy] *= -1.0;  // flip one spin and accept new spin config
                    MagneticMoment += (double) 2*SpinMatrix[ix][iy];
                    Energy += (double) deltaE;
                    AcceptedCycles++;
                }
            }
        }
        // update expectation values  for local node
        ExpectationValues[0] += Energy;
        ExpectationValues[1] += Energy*Energy;
        ExpectationValues[2] += MagneticMoment;
        ExpectationValues[3] += MagneticMoment*MagneticMoment;
        ExpectationValues[4] += fabs(MagneticMoment);
//        WriteResultstoFile(NSpins, cycles, Temperature, ExpectationValues, AcceptedCycles);
    }
} // end of Metropolis sampling over spins

// function to initialise energy, spin matrix and magnetization
void InitializeLattice(int NSpins, double** SpinMatrix,  double& Energy, double& MagneticMoment)
{
    // setup spin matrix and initial magnetization
    for(int x =0; x < NSpins; x++) {
        for (int y= 0; y < NSpins; y++){
            // Set RNG for lattice
            std::random_device rd;
            std::mt19937_64 gen(rd());
            std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
            if(RandomNumberGenerator(gen)<0.5){
                SpinMatrix[x][y] = 1.0; // spin orientation for the ground state
            } else {SpinMatrix[x][y] = -1.0;}
            MagneticMoment +=  (double) SpinMatrix[x][y];
        }
    }
    // setup initial energy
    for(int x =0; x < NSpins; x++) {
        for (int y= 0; y < NSpins; y++){
            Energy -=  (double) SpinMatrix[x][y]*
                    (SpinMatrix[PeriodicBoundary(x,NSpins,-1)][y] +
                    SpinMatrix[x][PeriodicBoundary(y,NSpins,-1)]);
        }
    }
}// end function initialise



void WriteResultstoFile(int NSpins, int MCcycles, double temperature, double* ExpectationValues, int AcceptedCycles)
{
    double norm = 1.0/((double) (MCcycles));  // divided by  number of cycles
    double E_ExpectationValues = ExpectationValues[0]*norm;
    double E2_ExpectationValues = ExpectationValues[1]*norm;
    double M_ExpectationValues = ExpectationValues[2]*norm;
    double M2_ExpectationValues = ExpectationValues[3]*norm;
    double Mabs_ExpectationValues = ExpectationValues[4]*norm;
    // all expectation values are per spin, divide by 1/NSpins/NSpins
    double Evariance = (E2_ExpectationValues- E_ExpectationValues*E_ExpectationValues)/NSpins/NSpins;
    double Mvariance = (M2_ExpectationValues - M_ExpectationValues*M_ExpectationValues)/NSpins/NSpins;
    double Mvarianceabs = (M2_ExpectationValues - Mabs_ExpectationValues*Mabs_ExpectationValues)/NSpins/NSpins;
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(15) << setprecision(8) << temperature;
    //ofile << setw(15) << setprecision(8) << MCcycles;
    //ofile << setw(15) << setprecision(8) << AcceptedCycles;
    ofile << setw(15) << setprecision(8) << E_ExpectationValues/NSpins/NSpins;
    ofile << setw(15) << setprecision(8) << Evariance/temperature/temperature;
    ofile << setw(15) << setprecision(8) << M_ExpectationValues/NSpins/NSpins;
    ofile << setw(15) << setprecision(8) << Mvariance/temperature;
    ofile << setw(15) << setprecision(8) << Mabs_ExpectationValues/NSpins/NSpins;
    ofile << setw(15) << setprecision(8) << Mvarianceabs/temperature << endl;
} // end output function

