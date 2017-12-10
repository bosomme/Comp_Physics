//
//  twodim.cpp
//  Project_5
//
//  Created by Børge Olsen-Hagen Sømme on 09.12.2017.
//  Copyright © 2017 Børge Olsen-Hagen Sømme. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#include "lib.h"

using namespace std;

ofstream ofile;

void initialize_matrix(int, double, double**);
void print_to_file_2D(string, double**, int, int, double);



int main(int argc, const char * argv[]) {
    string outfile_name;
    
    outfile_name = "2-D";
    
    int L = 1;
    int x_steps = 100; int y_steps = x_steps;
    double h = (double)L/x_steps; double dt = 0.5*h*h;
    
    
    // dx = 0.1 -> dt = 0.005
    // dx = 0.01 -> dt = 0.00005
    
    
    int time_steps = 10*(floor(1/(dt*10)+0.1));
    
    dt = 1e-5; time_steps = 1/dt+1; cout << time_steps;
    
    double alpha = dt/(h*h); //alpha = 0.5?
    double minus4alpha = 1-4*alpha;
    
    
    double **u; u = (double**) matrix(x_steps+1, y_steps+1, sizeof(double));
    double **u_new; u_new = (double**) matrix(x_steps+1, y_steps+1, sizeof(double));
    
    initialize_matrix(x_steps+1, h, u); initialize_matrix(x_steps+1, h, u_new);
    
    print_to_file_2D(outfile_name, u, x_steps, y_steps, (double)0/time_steps);
    
    for(int k=1; k<time_steps; k++){
        // twodimensional(x_steps+1, y_steps+1, alpha, u, u_new);
        
        // u[i, j, k+1] = u[i, j, k] + alpha*(u[i+1,j,k] + u[i-1,j,k] + u[i,j+1,k] + u[i,j-1,k] - 4u[i,j,k])
        for(int i=1; i<x_steps; i++){
            for(int j=1; j<y_steps; j++){
                // u_new[i][j] = u[i+1][j] + alpha*(u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1] - 4*u[i][j]);
                u_new[i][j] = minus4alpha*u[i][j] + alpha*(u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1]);
            }
        }
        u_new[0][0] = u_new[0][y_steps] = u_new[x_steps][0] = u_new[x_steps][y_steps] = 0;
        for(int i=0; i<x_steps+1; i++){
            for(int j=0; j<y_steps+1; j++){
                u[i][j] = u_new[i][j];
            }
        }
        
        
        // Print to file
        if (k==0.05*time_steps){
            print_to_file_2D(outfile_name, u, x_steps, y_steps, (double)k/time_steps);
        }
        else if (k==0.5*time_steps){
            print_to_file_2D(outfile_name, u, x_steps, y_steps, (double)k/time_steps);
        }

    }
    
    // Freeing memory
    free_matrix((void **) u);
    free_matrix((void **) u_new);
    
    return 0;
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
