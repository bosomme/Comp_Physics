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
#include <string>

#include "lib.h"


using namespace std;

ofstream ofile;


void print_to_file(string, double*, int, double);

void Forward_Euler(int, double, double*, double*);
void Backward_Euler(int, double, double*, double*);
void Crank_Nicolson(int, double, double*);

void twodimensional(int, int, double, double**, double**);
void initialize_matrix(int, double, double**);
void print_to_file_2D(string, double**, int, int, double);


int main(int argc, const char * argv[]) {
    string outfile_name;
    
    int L = 1;
    int x_steps = 100; int y_steps = x_steps;
    double dx = (double)L/x_steps; double dt = 0.5*dx*dx;
    
    // dx = 0.1 -> dt = 0.005
    // dx = 0.01 -> dt = 0.00005
    
    int time_steps = 10*(floor(1/(dt*10)+0.1));
    double alpha = dt/(dx*dx); //alpha = 0.5?
    
    double *u = new double[x_steps+1];
    double *u_new = new double[x_steps+1];
    
    
    // Decide method
    cout << "Please give wanted method number:" << endl;
    cout << "1 = Forward Euler, 2 = Backward Euler, 3 = Crank-Nicolson, 4 = Two-dimensional" << endl;
    int method;
    cin >> method;
    
    for(int i=0; i<=x_steps; i++) u[i] = u_new[i] = 0;
    u[x_steps] = u_new[x_steps] = 1;
    
    
    for(int j=1; j<time_steps; j++){
        if (method == 1){
            Forward_Euler(x_steps+1, alpha, u, u_new);
            outfile_name = "F_E";
        }
        else if (method == 2){
            Backward_Euler(x_steps+1, alpha, u, u_new);
            outfile_name = "B_E";
        }
        else if (method == 3){
            Crank_Nicolson(x_steps+1, alpha, u);
            outfile_name = "C-N";
        }
        else if (method == 4){
            // delete [] u; delete [] u_new;
            double **u; u = (double**) matrix(x_steps+1, y_steps+1, sizeof(double));
            double **u_new; u_new = (double**) matrix(x_steps+1, y_steps+1, sizeof(double));
            
            initialize_matrix(x_steps+1, dx, u); initialize_matrix(x_steps+1, dx, u_new);
            
            twodimensional(x_steps+1, y_steps+1, alpha, u, u_new);
            outfile_name = "2-D";
            
            // Print to file (Twodimensional)
            if (j==0.05*time_steps){
                print_to_file_2D(outfile_name, u, x_steps, y_steps, (double)j/time_steps);
            }
            else if (j==0.5*time_steps){
                print_to_file_2D(outfile_name, u, x_steps, y_steps, (double)j/time_steps);
            }
        }
        
        
        // Print to file (Onedimensional)
        if (method==1 || method==2 || method==3){
            if (j==0.05*time_steps){
                print_to_file(outfile_name, u, x_steps, (double)j/time_steps);
            }
            else if (j==0.5*time_steps){
                print_to_file(outfile_name, u, x_steps, (double)j/time_steps);
            }
        }
    }
    
    // Freeing memory
    if (method==1 || method==2 || method==3) {delete [] u; delete [] u_new;}
    else if (method==4) {free_matrix((void **) u); free_matrix((void **) u_new);}
    
    return 0;
}

void twodimensional(int Number_of_XSteps, int Number_of_YSteps, double alpha, double** u, double** u_new){
    // u[i, j, k+1] = u[i, j, k] + alpha*(u[i+1,j,k] + u[i-1,j,k] + u[i,j+1,k] + u[i,j-1,k] - 4u[i,j,k])
    for(int i=1; i<Number_of_XSteps-1; i++){
        for(int j=1; j<Number_of_YSteps-1; j++){
            u_new[i][j] = u[i+1][j] + alpha*(u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1] - 4*u[i][j]);
        }
    }
    u_new[0][0] = u_new[0][Number_of_XSteps-1] = u_new[Number_of_XSteps-1][0] = u_new[Number_of_XSteps-1][Number_of_XSteps-1] = 0;
    for(int i=1; i<Number_of_XSteps; i++){
        for(int j=1; j<Number_of_YSteps; j++){
            u[i][j] = u_new[i][j];
        }
    }
}


void initialize_matrix(int n, double dx, double **A){
    double pi = M_PI;
    double x; double y;
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            if(i == 0 || i == n-1 || j == 0 || j == n-1){
                A[i][j] = 0;
            }
            else{
                x = i*dx; y = j*dx;
                A[i][j] = sin(pi*x)*sin(pi*y);
            }
        }
    }
}

void print_to_file_2D(string outfile_name, double **u, int nx, int ny, double t){
    string outfile = outfile_name;
    outfile.append("_t="); outfile.append(to_string(t));
    outfile.append("_dx="); outfile.append(to_string(1.0/nx));
    outfile.append(".txt");
    ofile.open(outfile);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    for(int i=0; i<=nx; i++){
        ofile << (double) i/nx << "  ";
        for(int j=0; j<=ny; j++){
            ofile << setprecision(8) << u[i][j] << "  ";
        }
        ofile << endl;
    }
    ofile.close();
}






void print_to_file(string outfile_name, double *u, int n, double t){
    string outfile = outfile_name;
    outfile.append("_t="); outfile.append(to_string(t));
    outfile.append("_dx="); outfile.append(to_string(1.0/n));
    outfile.append(".txt");
    ofile.open(outfile);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    for(int i=0; i<=n; i++){
        ofile << (double) i/n << "  ";
        ofile << setprecision(8) << u[i] << endl;
    }
    ofile.close();
}


// Function that solves Au_new = u, with A being a tridiagonal matrix.
void tridiag_solver(int n, double *a, double *b, double *c, double *u, double *u_new){
    double *b_tilde = new double[n];
    double *u_tilde = new double[n];
    
    b_tilde[0] = b[0];
    u_tilde[0] = u[0];
    
    double temp;
    // Forward substitution
    for (int i=1; i<n; i++){
        temp = a[i-1]/b_tilde[i-1];
        b_tilde[i] = b[i] - c[i-1]*temp;
        u_tilde[i] = u[i] - u_tilde[i-1]*temp;
    }
    
    u_new[n-1] = u_tilde[n-1]/b_tilde[n-1];
    
    // Backward substitution
    for (int i=n-2; i>=0; i--){
        u_new[i] = (u_tilde[i] - c[i]*u_new[i+1])/b_tilde[i];
    }
    
    delete [] b_tilde;
    delete [] u_tilde;
}


// Forward Euler Algorithm:
// Function for computing one time step of PDE using the Forward Euler algorithm.
void Forward_Euler(int Number_of_XSteps, double alpha, double *u, double *u_new){
    // Compute solution for next time step
    for(int i=1; i<Number_of_XSteps; i++){
        u_new[i] = alpha*u[i-1] + (1-2*alpha)*u[i] + alpha*u[i+1];
    }
    
    // Boundary conditions
    u_new[0] = 0; u_new[Number_of_XSteps-1] = 1;
    
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
    for(int i=0; i<Number_of_XSteps; i++){
        a[i] = -alpha;      // Subdiagonal (same on each side)
        b[i] = 1 + 2*alpha; // Diagonal
    }
    
    u[Number_of_XSteps-1] += alpha;
    
    // Compute next time step
    tridiag_solver(Number_of_XSteps, a, b, a, u, u_new);
    
    // Boundary conditions
    u_new[0] = 0;
    u_new[Number_of_XSteps-1] = 1;
    
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
    utemp[0] = 0; utemp[Number_of_XSteps-1] += 2*alpha;
    
    // Solving equation:
    double* a = new double[Number_of_XSteps];
    double* b = new double[Number_of_XSteps];
    // Fill arrays for diagonal and subdiagonal elements:
    for(int i=0; i<Number_of_XSteps; i++){
        a[i] = -alpha;      // Subdiagonal (same on each side)
        b[i] = 2 + 2*alpha; // Diagonal
    }
    
    // Compute next time step:
    tridiag_solver(Number_of_XSteps, a, b, a, utemp, u);
    
    // Boundary conditions
    u[0] = 0;
    u[Number_of_XSteps-1] = 1;
    
    
    delete [] utemp; delete [] a; delete [] b;
}





