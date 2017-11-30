//
//  main.cpp
//  Project_5
//
//  Created by Børge Olsen-Hagen Sømme on 30.11.2017.
//  Copyright © 2017 Børge Olsen-Hagen Sømme. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <iomanip>


#include "lib.h"



void Forward_Euler(int, double, double*, double*);
void Backward_Euler(int, double, double*, double*);
void Crank_Nicolson(int, double, double*);

int main(int argc, const char * argv[]) {
    int L = 1;
    int x_steps = 100;
    int time_steps = 100;
    double dx = (double)L/x_steps; double dt = 0.5*dx*dx;
    
    double alpha = dt/(dx*dx); //alpha = 0.5?
    
    double *u = new double[x_steps];
    double *u_new = new double[x_steps];
    
    
    for(int j=1; j<time_steps; j++){
        // Forward_Euler(x_steps, alpha, u, u_new);
        // Backward_Euler(x_steps, alpha, u, u_new);
        Crank_Nicolson(x_steps, alpha, u);
    }
    
    
    
    delete [] u;
    delete [] u_new;
    
    return 0;
}

void tridiag_solver(int n, double *a, double *b, double *c, double *u, double *u_new){
    double *b_tilde = new double[n];
    double *u_tilde = new double[n];
    
    b_tilde[0] = b[0];
    u_tilde[0] = u[0];
    
    double temp = 0.0;
    
    // Forward substitution
    for (int i=1; i<n; i++)
    {
        temp = a[i-1]/b_tilde[i-1];
        b_tilde[i] = b[i] - c[i-1]*temp;
        u_tilde[i] = u[i] - u_tilde[i-1]*temp;
    }
    
    u_new[n-1] = u_tilde[n-1]/b_tilde[n-1];
    
    // Backward substitution
    for (int i=n-2; i>=0; i--)
    {
        u_new[i] = (u_tilde[i] - c[i]*u_new[i+1])/b_tilde[i];
    }
    
    
    delete [] b_tilde; delete [] u_tilde;
}

// Forward Euler Algorithm:
// Function for computing one time step of PDE using the Forward Euler algorithm.
void Forward_Euler(int Number_of_XSteps, double alpha, double *u, double *u_new){
    // Compute solution for next time step
    for(int i=1; i<Number_of_XSteps; i++){
        u_new[i] = alpha*u[i-1] + (1-2*alpha)*u[i] + alpha*u[i+1];
    }
    // Update time solution with the new one
    for(int i=1; i<Number_of_XSteps; i++){
        u[i] = u_new[i];
    }
}

// Backward Euler Algorithm:
//
void Backward_Euler(int Number_of_XSteps, double alpha, double *u, double *u_new){
    double* a = new double[Number_of_XSteps];
    double* b = new double[Number_of_XSteps];
    // Fill arrays for diagonal and subdiagonal elements:
    for(int i=1; i<Number_of_XSteps; i++){
        a[i] = -alpha;      // Subdiagonal (same on each side
        b[i] = 1 + 2*alpha; // Diagonal
    }
    
    // Compute next time step:
    tridiag_solver(Number_of_XSteps, a, b, a, u, u_new);
    
    // Update time solution with the new one
    for(int i=1; i<Number_of_XSteps; i++){
        u[i] = u_new[i];
    }
    
    
    delete [] a; delete [] b;
}


// Crank-Nicolson Scheme:
//
void Crank_Nicolson(int Number_of_XSteps, double alpha, double *u){
    // Right-hand side:
    double* utemp = new double[Number_of_XSteps];
    for(int i=1; i<Number_of_XSteps; i++){
        utemp[i] = alpha*u[i-1] + (2-2*alpha)*u[i] + alpha*u[i+1];
    }
    
    // Solving equation:
    double* a = new double[Number_of_XSteps];
    double* b = new double[Number_of_XSteps];
    // Fill arrays for diagonal and subdiagonal elements:
    for(int i=1; i<Number_of_XSteps; i++){
        a[i] = -alpha;      // Subdiagonal (same on each side
        b[i] = 2 + 2*alpha; // Diagonal
    }
    
    // Compute next time step:
    tridiag_solver(Number_of_XSteps, a, b, a, utemp, u);
    
    
    
    delete [] utemp; delete [] a; delete [] b;
}





