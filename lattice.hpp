//
//  lattice.hpp
//  exam
//
//  Created by Solfrid Johansen on 04/05/2020.
//  Copyright © 2020 Solfrid Johansen. All rights reserved.
//

#pragma once

// Increases matrix size allowed to allocate (might not be needed, depends on system)
#define ARMA_64BIT_WORD 1

#include <iostream>
#include <armadillo>
#include <string>
#include <vector>
#include <cmath>
#include <complex>

using namespace std;
using namespace arma;

class Lattice{
    private:
        int N_x; // Size of lattice in x-direction
        int N_y; // Size of lattice in y-direction
        Mat<double> latticeOne;
        int J;
        int k_b;
        double n1;
        double T; // [K] temperature
        int nr_sweeps; // Nr. sweeps after equilibrium is reached
        int equilibrium; // Nr. sweeps until equilibrium is reached
        vec conv_test;
    public:

    Lattice(int N_x, int N_y){
        this->N_x = N_x;
        this->N_y = N_y;
        mon_jasnow();
        J = 1;
        k_b = 1;
        // Må ha lav temp, enda lavere for å se ordering i plot
        T = 1.9;
        equilibrium = 1000;
        nr_sweeps = 1000;
        n1 = 1.0/(nr_sweeps*N_x*N_x);
    }

    void mon_jasnow(){
        // Create lattice of correct size, and fill with random spins (either -1 or 1)
        latticeOne.set_size(N_y,N_x+2);
        
        // Fills with random values
        latticeOne = round(randu(N_y, N_x + 2));
        // Values less than 1 is set to -1
        latticeOne(find(latticeOne<1)).fill(-1);
        // Fixed spin boundary conditions
        latticeOne.col(0).fill(1);
        latticeOne.col(N_x+1).fill(1);

    }

    void extended_mon_jasnow();
    
    double getMagnetization();
    
    void mcSweep(int method);
    
    int mod(int a, int b);
    
    double computedeltaH(int i, int j);
    double computedeltaHExtended(int i, int j);
    
    double computeTorusHamiltonian();
    double computeKleinHamiltonian();
    
    bool tryFlip(int i, int j, int method);
    
    double metropolis();
    double metropolis_extended();
    
    // Functions used to get data needed for plotting
    void tauOfTemp();
    double tau_o_plot();
    void convergence();
    void computeNtau();
};
