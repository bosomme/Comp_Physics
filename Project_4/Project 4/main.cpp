//
//  main.cpp
//  Project 4
//
//  Created by Børge Olsen-Hagen Sømme on 06.11.2017.
//  Copyright © 2017 Børge Olsen-Hagen Sømme. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include "lib.h"

using namespace std;

ofstream ofile;

void initialize(int, int**, double&, double&, double, int);
void Metropolis_algorithm(long&, int, int**, double&, double&, double*);


inline int periodic(int i, int limit, int add){
    return (i+limit+add) % (limit);
}


int main(int argc, const char * argv[]) {
    // J = 1, k_B = 1
    
    ofile.open("L=80_1E6.txt");
    
    ofile << "Energies and Magnetization as function of temperature" << endl;
    ofile << "Temperature - E - M - C_V - X" << endl;
    
    /*
    ofile << "Mean energy and magnetization as function of number of MonteCarlo cycles" << endl;
    ofile << "Number of MC - E - M" << endl;
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    */
    int random = 1; // If equal to 1 set up a random starting matrix, if other value set up initial matrix with all spins up
    
    
    long idum = -1;
    // Initialize matrix:
    int number_of_spins = 80; // Number of spins (matrix-points) in each direction (L)
    
    int **spin_lattice;
    double w[17], average[5], E, M;
    // average[0] = E, average[1] = E^2, average[2] = M, average[3] = |M|^2, average[4] = |M|
    
    double initial_temperature = 2.1;
    double final_temperature = 2.3;
    double temp_stepsize = 0.05;
    
    int number_of_MC = 1E6;
    

    double heat_capacity, susceptibility;
    double E_out, EE_out, M_out, MM_out, absM_out;
    
    
    double norm = 1/(double) number_of_MC;
    
    
    spin_lattice = (int**) matrix(number_of_spins, number_of_spins, sizeof(int));
    
    
    for(double T=initial_temperature; T<=final_temperature; T+=temp_stepsize){
        E = M = 0.;
        initialize(number_of_spins, spin_lattice, E, M, T, random);
        
        for(int dE=-8; dE<=8; dE++){
            w[dE+8] = 0;
        }
        for(int dE=-8; dE<=8; dE+=4){
            w[dE+8] = exp(-dE/T);
        }
        for(int i=0; i<5; i++){
            average[i] = 0.;
        }
        
        for(int i=1; i<= number_of_MC; i++){
            Metropolis_algorithm(idum, number_of_spins, spin_lattice, E, M, w);
            
            average[0] += E; average[1] += E*E;
            average[2] += M; average[3] += M*M; average[4] += abs(M);
            
            
            norm = 1/(double)i;
        }
        
        
        E_out = average[0]*norm; EE_out = average[1]*norm;
        M_out = average[2]*norm; MM_out = average[3]*norm; absM_out = average[4]*norm;
        
        heat_capacity = (EE_out - E_out*E_out)/(number_of_spins*number_of_spins);
        susceptibility = (MM_out - absM_out*absM_out)/(number_of_spins*number_of_spins);
        ofile << T << "   " << E_out << "   " << absM_out << "   ";
        ofile << heat_capacity << "   " << susceptibility  << endl;
        
        
    }
    
    ofile.close();
    free_matrix((void**) spin_lattice);
}


void initialize(int number_of_spins, int **spin_lattice, double &E, double &M, double T, int random){
    if (random == 1) {
        long idum_i = -1;
        for(int x=0; x<number_of_spins; x++){
            for(int y=0; y<number_of_spins; y++){
                spin_lattice[y][x] = 1;
                if (ran1(&idum_i) < 0.5) spin_lattice[y][x] *= -1;
                M += (double) spin_lattice[y][x];
            }
        }
    }
    else {
        for(int x=0; x<number_of_spins; x++){
            for(int y=0; y<number_of_spins; y++){
                spin_lattice[y][x] = 1;
                M += (double) spin_lattice[y][x];
            }
        }
    }
    for(int x=0; x<number_of_spins; x++){
        for(int y=0; y<number_of_spins; y++){
            E -= (double) spin_lattice[y][x]*
                 (spin_lattice[y][periodic(x, number_of_spins, -1)]
                 + spin_lattice[periodic(y, number_of_spins, -1)][x]);
        }
    }
}


void Metropolis_algorithm(long &idum, int number_of_spins, int **spin_lattice, double &E, double &M, double *w) {
    for(int x=0; x<number_of_spins; x++){
        for(int y=0; y<number_of_spins; y++){
            int ran_x = (int) (ran1(&idum)*(double)number_of_spins);
            int ran_y = (int) (ran1(&idum)*(double)number_of_spins);
            
            int delta_E = 2*spin_lattice[ran_y][ran_x]*
            (spin_lattice[ran_y][periodic(ran_x, number_of_spins, -1)] +
             spin_lattice[periodic(ran_y, number_of_spins, -1)][ran_x] +
             spin_lattice[ran_y][periodic(ran_x, number_of_spins, 1)] +
             spin_lattice[periodic(ran_y, number_of_spins, 1)][ran_x]);
            
            if (ran1(&idum) <= w[delta_E+8]){
                spin_lattice[ran_y][ran_x] *= -1;
                M += (double) 2*spin_lattice[ran_y][ran_x];
                E += (double) delta_E;
            }
        }
    }
}


