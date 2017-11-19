#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <string>
#include "ising.h"
#include "mpi.h"
using namespace  std;
// output file
ofstream ofile;

int main(int argc, char* argv[])
{
    string filename = "IsingTempRange";
    int NSpins=100;
    int MCcycles=1e6; //Maximum number of MC cycles
    int CutOff = 5e4;
    double InitialTemp=2.0;   //kT/J
    double Temperature = InitialTemp;
    double FinalTemp=2.3;
    double TempStep=0.005;
    int numprocs, my_rank;

    MPI_Init(&argc,&argv);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

 if(my_rank==0){
        string fileout = filename;
        string argument = to_string(NSpins);
        fileout.append(argument);
        fileout.append(".txt");
        ofile.open(fileout);
   }

    for (Temperature = InitialTemp; Temperature <= FinalTemp; Temperature+=TempStep){
        double* ExpectationValues = new double [6];
        for(int ii=0; ii<=5; ii++) ExpectationValues[ii]=0.;

        MetropolisSampling(NSpins, MCcycles, Temperature, ExpectationValues, CutOff);

        double* TotalExpectationValues = new double [6];
        for(int ii=0; ii<=5; ii++){
            TotalExpectationValues[ii]=0.;
            MPI_Reduce(&ExpectationValues[ii], &TotalExpectationValues[ii], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }
        WriteResultstoFile(NSpins, (MCcycles-CutOff)*numprocs, Temperature, TotalExpectationValues, 0);
        delete [] ExpectationValues;
       // delete [] TotalExpectationValues;
    }
    if(my_rank==0)
        ofile.close();  // close output file

    MPI_Finalize();
    return 0;
}



