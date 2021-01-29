//
//  lattice.cpp
//  exam
//
//  Created by Solfrid Johansen on 04/05/2020.
//  Copyright © 2020 Solfrid Johansen. All rights reserved.
//

#include "lattice.hpp"

// Computes modulo (not just remainder, like the % operator does)
int Lattice::mod(int a, int b)
{
    int r = a % b;
    return r < 0 ? r + b : r;
}

double Lattice::computedeltaH(int i, int j){
    // Periodic boundary  conditions in y-direction (top and bottom rows)
    // Modulo restricts values between 0 and N_y-1
    int bottom = latticeOne(mod(i+1, N_y), j);
    int top  = latticeOne(mod(i-1, N_y), j);
    int left = latticeOne(i, j-1);
    int right = latticeOne(i, j+1);
    int nearestNeighbors = top + bottom + left + right;
    double deltaH = 2*latticeOne(i,j)*nearestNeighbors;
    return deltaH;
}

double Lattice::computedeltaHExtended(int i, int j){
    // Periodic boundary  conditions in y-direction (top and bottom rows)
    // Modulo restricts values between 0 and N_y-1
    int bottom = latticeOne(mod(i+1, N_y), j);
    int top  = latticeOne(mod(i-1, N_y), j);
    int left = latticeOne(i, mod(j-1, N_x));
    int right = latticeOne(i, mod(j+1,N_x));
    int nearestNeighbors = top + bottom + left + right;
    double deltaH = 2*latticeOne(i,j)*nearestNeighbors;
    return deltaH;
}


bool Lattice::tryFlip(int i, int j, int method){
    double deltaH = 0;
    if(method==0){
        deltaH = computedeltaH(i, j);
    }else{
        deltaH = computedeltaHExtended(i, j);
    }

    // Random integer between 0 and 1
    double random = randu();
    bool flip = false;
    if(deltaH < 0){
        flip = true;
    }else{
        double w = exp(-deltaH/T);
        if(w > random){
            flip = true;
        }
    }
    return flip;

}

double Lattice::getMagnetization(){
    // Magnetization
    return accu(latticeOne);
}

void Lattice::mcSweep(int method){
    // For each sweep, try to flip N^2 times
    for(int k = 0; k < N_x*N_x; k++){
        double delta_m = 0;
        int random_i = rand()%N_y;
        int random_j = 0;
        if(method==0){
            random_j = rand()%(N_x) + 1;
        }
        else{
            random_j = rand()%N_x;
        }
        bool flip = tryFlip(random_i, random_j, method);
        if(flip){
            // To flip spin, multiply by minus one.
            latticeOne(random_i, random_j) *= -1;
        }
    }
}

double Lattice::metropolis(){
    // Compute magnetization for this sweep
    double sweep_m = 0;
    // To perform convergence test
    //conv_test.set_size(equilibrium + nr_sweeps);
    // Must reach equilibrium
    for(int t = 0; t < equilibrium; t++){
        mcSweep(0);
        // UNCOMMENT FOR CONVERGENCE TEST
        //double mag = getMagnetization();
        //conv_test(t) = mag;
    }
    for(int t = 0; t < nr_sweeps; t++){
        mcSweep(0);
        double mag = getMagnetization();
        sweep_m += mag;
    }
    // Returns the magnetization
    return sweep_m*n1;
}

double Lattice::metropolis_extended(){
    // Compute magnetization for this sweep
    double H_t = 0;
    double H_k = 0;
    double sweep_m = 0;

    // UNCOMMENT TO SAVE STATE OF LATTICE BEFORE SWEEP
    //latticeOne.save("../data/lattice_before.csv", csv_ascii);

    // Must reach equilibrium
    for(int t = 0; t < equilibrium; t++){
        mcSweep(1);

    }
    for(int t = 0; t < nr_sweeps; t++){
        mcSweep(1);
        H_t = computeTorusHamiltonian();
        H_k = computeKleinHamiltonian();
        double m = (H_k - H_t);
        sweep_m += exp(-m / T);
    }
    // UNCOMMENT TO SAVE STATE OF LATTICE AFTER SWEEP
    //latticeOne.save("../data/lattice_after.csv", csv_ascii);
    return sweep_m / nr_sweeps;
}

void Lattice::tauOfTemp(){
    // compute up to T_c + 1
    double T_c = 2 / log(1 + sqrt(2)) + 1;
    // Start temperature
    T = 0.1;
    int nr_computations = 15;

    double step = (T_c - T) / nr_computations;

    vec temp(nr_computations);
    vec taus(nr_computations);

    int count = 0;
    while(count < nr_computations){
        // computes the magnetization
        
        // USE BELOW FOR ORDINARY
        //mon_jasnow();
        //taus(count) = 2*(metropolis());

        // USE BELOW FOR EXTENDED
        extended_mon_jasnow();
        taus(count) = - T / N_y * log(metropolis_extended());
    
        temp(count) = T;
        cout << count << endl;
        count++;
        T+=step;
    }
    // CHANGE FILENAMES APPROPRIATELY
    temp.save("../data/temps_extended_100.csv", csv_ascii);
    taus.save("../data/tau_extended_100.csv", csv_ascii);
}

void Lattice::computeNtau(){
  
    double T_c = 2/log(1+sqrt(2));
    // Temperature
    T = T_c;

    int N_min = 5;
    int N_max = 25;
    vec tauN(N_max-N_min);
    vec N(N_max-N_min);
    for(int i = N_min; i < N_max; i++){
        N_x = i;
        N_y = i;
        n1 = 1.0/(nr_sweeps*N_x*N_x);
        
        //----------------------------
        // FOR ORDINARY
        //mon_jasnow();
        //tauN(i-N_min) = 2*metropolis()*i;
        //----------------------------

        //---------------------------
        // FOR EXTENDED
        double tau = 0;
        for(int j = 0; j < 100; j++){
            arma_rng::set_seed_random();
            extended_mon_jasnow();
            tau +=  - T / N_y * log(metropolis_extended());
        }
        //--------------------------
        tau = tau/10;
        tauN(i-N_min) = tau*i;
        N(i-N_min) = i;
        cout << i << endl;
    }
    tauN.save("../data/tauN_ext.csv", csv_ascii);
    N.save("../data/N_ext.csv", csv_ascii);
}

void Lattice::convergence(){
    // NB: UNCOMMENT SECTIONS IN metropolis() BEFORE
    T = 1.9;
    mon_jasnow();
    metropolis();
    conv_test.save("../data/conv_200.csv", csv_ascii);
}

double Lattice::computeTorusHamiltonian(){
    double H_t = 0;
    for(int i = 0; i < N_x; i++) {
        for(int j = 0; j < N_y; j++) {
            // Torus
            int bottom = latticeOne(mod(i + 1, N_y), j);
            int top = latticeOne(mod(i - 1, N_y), j);
            int left = latticeOne(i, mod(j - 1, N_x));
            int right = latticeOne(i, mod(j + 1, N_x));
            int nearestNeighbors = top + bottom + left + right;
            H_t += -latticeOne(i, j) * nearestNeighbors;
        }
    }
    // Avoid overcounting
    return H_t*0.5;
}

double Lattice::computeKleinHamiltonian(){
    double H_k = 0;
    for(int i = 0; i < N_x; i++) {
        for(int j = 0; j < N_y; j++) {
            int point = latticeOne(i,j);
            // Klein bottle structure of lattice
            int bottom = latticeOne(mod(i + 1, N_y), j);
            int top = latticeOne(mod(i - 1, N_y), j);
            int left = 0;
            int right = 0;
            // At leftmost column
            if(j == 0){
                left = -latticeOne(N_x-1-i, N_y-1);
                right  = latticeOne(i, j+1);

            } // At rightmost column
            else if(j == N_x-1){
                right = -latticeOne(N_x-1-i, 0);
                left = latticeOne(i, j-1);
            }else{
                left = latticeOne(i, j-1);
                right = latticeOne(i, j+1);
            }

            int nearestNeighbors = top + bottom + left + right;
            H_k += -latticeOne(i, j) * nearestNeighbors;
        }
    }
    return H_k*0.5;
}

void Lattice::extended_mon_jasnow(){
    // Set size and fill with random spins (either up or down)
    latticeOne.set_size(N_y,N_x);
    latticeOne = round(randu(N_y, N_x));
    latticeOne(find(latticeOne<1)).fill(-1);
}

double Lattice::tau_o_plot(){
    // Select the temperature to plot for, between 0.1 and 2.23 (T_c)
    T = 2.23; // [K]

    double T_c = 2/log(1+sqrt(2));
    double t = (T_c-T)/T_c;

    int N_min = 2;
    int N_max = 150;

    vec tau_temp(N_max-N_min);
    vec inv_N_t(N_max-N_min);

    // Go thorugh for different T
    for(int i = N_min; i < N_max; i++){
        N_x = i;
        N_y = i;
        n1 = 1.0/(nr_sweeps*N_x*N_x);
        double tau = 0;
        double nr_it = 1;
        for(int k= 0; k < nr_it; k++){
            // Compute nr_it times per lattice size
            extended_mon_jasnow();
            tau += - T / N_y * log(metropolis_extended());
        }
        tau = tau/nr_it;
        tau_temp(i-N_min) = tau;
        cout <<  tau_temp(i-N_min) << endl;
        inv_N_t(i-N_min) = i;
        cout << i << endl;

    }
    cout << t << endl;
    tau_temp.save("../data/tau_223.csv", csv_ascii);
    inv_N_t.save("../data/N_223.csv", csv_ascii);
}


