The code is divided into two files:
- main.cpp: contains the main functions
- ising.cpp: contains all the functions called in main

Set initial parameters in the first lines of "main.h"
Comment/uncomment printing statement depending on whether you want to get data after very MC cycle or after the number of cycles in the variable "MCcycles" (for example when running through the temperature interval comment the statement in "ising.cpp" and uncomment the one in "main.cpp")
The output file is organized in data columns: temperature | MC cycles | Accepted Cycles | MeanEnergy | specific heat MeanMagnetization | susceptibility | Mean Absolute Magnetization | susceptibility(computed with mean absolute magnetization) 
Note that all the data are normalized by the total number of spins.
If the function "PrintEnergyCounts" is uncommented, it will print on another file the counts of possible energy values to analyze their distribution.
 
The parallelized version of the program has been compiled in the machines at the university lab with the following commands:
- Compile: /usr/lib64/openmpi/bin/mpic++ -std=c++11 -o3 parallel_main.cpp -o ising
- Run: /usr/lib64/openmpi/bin/mpirun -np 4 ising
