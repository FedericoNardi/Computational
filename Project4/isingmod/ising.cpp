#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <string>
#include "ising.h"
using namespace  std;
extern ofstream ofile;

void MetropolisSampling(int NSpins, int MCcycles, double Temperature, double* ExpectationValues)
{
    // Initialize RNG
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
    // Initialize the lattice spin values
    double** SpinMatrix = new double*[NSpins];
    for(int i = 0; i<NSpins; i++){
        SpinMatrix[i]=new double[NSpins];
    }
    double Energy = 0.;     double MagneticMoment = 0.;
    InitializeLattice(NSpins, SpinMatrix, Energy, MagneticMoment);

    double* EnergyDifference = new double[5];
    for( int de =-8; de <= 8; de+=4){
        double index = 0.25*(de+8);
        EnergyDifference[(int)index] = exp(-de/Temperature);
    }


    int AcceptedMoves=0;
    int EnergyCounts[401];
        for(int i=0; i<401; i++) EnergyCounts[i]=0;
    // Monte Carlo cycles start
    for (int cycles = 1; cycles <= MCcycles; cycles++){
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
                    AcceptedMoves++;
                    EnergyCounts[(int)((Energy+800)/4)]++;
                }
            }
        }
        // update (local) expectation values
//        if(cycles>CutOff){
            Update(ExpectationValues, Energy, MagneticMoment);
            WriteResultstoFile(NSpins, cycles, Temperature, ExpectationValues, AcceptedMoves);
//        }
    }
//    // Normalize Energy counts to get probabilities
//    for(int k=0; k<401; k++) EnergyCounts[k] = EnergyCounts[k]/MCcycles;
//    // Print Energy counts on file
//    PrintEnergyCounts(EnergyCounts, Temperature);
//    WriteResultstoFile(NSpins, MCcycles, Temperature, ExpectationValues, AcceptedMoves);
}

void PrintEnergyCounts(int EnergyCounts[], double Temperature){
    ofstream CountsFile;
    string ofname = "EnergyCounts_T";
    ofname.append(to_string((int)Temperature));
    ofname.append(".txt");
    CountsFile.open(ofname);
    CountsFile << setiosflags(ios::showpoint | ios::uppercase);
    for (int k=0; k<401; k++) CountsFile << setw(15) << setprecision(8) <<EnergyCounts[k] <<"\n";
    CountsFile.close();
}

void Update(double* ExpectationValues, double Energy, double MagneticMoment){
    ExpectationValues[0] += Energy;
    ExpectationValues[1] += Energy*Energy;
    ExpectationValues[2] += MagneticMoment;
    ExpectationValues[3] += MagneticMoment*MagneticMoment;
    ExpectationValues[4] += fabs(MagneticMoment);
}


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
            } else SpinMatrix[x][y] = -1.0;
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
}



void WriteResultstoFile(int NSpins, int MCcycles, double temperature, double* ExpectationValues, int AcceptedCycles)
{
    double norm = 1.0/((double) (MCcycles));  // divided by  number of cycles
    double MeanEnergy = ExpectationValues[0]*norm;
    double MeanEnergySquared = ExpectationValues[1]*norm;
    double MeanMagMoment = ExpectationValues[2]*norm;
    double MeanMagMomentSquared = ExpectationValues[3]*norm;
    double MeanAbsMagMoment = ExpectationValues[4]*norm;
    // all expectation values are per spin, divide by 1/NSpins/NSpins
    double Evariance = (MeanEnergySquared- MeanEnergy*MeanEnergy)/NSpins/NSpins;
    double Mvariance = (MeanMagMomentSquared - MeanMagMoment*MeanMagMoment)/NSpins/NSpins;
    double Mvarianceabs = (MeanMagMomentSquared - MeanAbsMagMoment*MeanAbsMagMoment)/NSpins/NSpins;
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(15) << setprecision(8) << temperature;
    ofile << setw(15) << setprecision(8) << MCcycles;
    ofile << setw(15) << setprecision(8) << AcceptedCycles;
    ofile << setw(15) << setprecision(8) << MeanEnergy/NSpins/NSpins;
    ofile << setw(15) << setprecision(8) << Evariance/temperature/temperature;
    ofile << setw(15) << setprecision(8) << MeanMagMoment/NSpins/NSpins;
    ofile << setw(15) << setprecision(8) << Mvariance/temperature;
    ofile << setw(15) << setprecision(8) << MeanAbsMagMoment/NSpins/NSpins;
    ofile << setw(15) << setprecision(8) << Mvarianceabs/temperature << "\n";
}

