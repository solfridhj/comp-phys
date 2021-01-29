//

//
//  Created by Solfrid Johansen on 04/05/2020.
//  Copyright © 2020 Solfrid Johansen. All rights reserved.
//

#include "lattice.hpp"
#include <iostream>

int main(int argc, char *argv[]) {
    cout << "Welcome to Mon-Jasnow! Go to main.cpp and uncomment to run programs." << endl;
    // Set seed
    arma_rng::set_seed(10);
    Lattice lat(100, 100);
    //lat.extended_mon_jasnow();
    //lat.metropolis_extended();
    //lat.computeTau();
    //lat.tauOfTemp();
    //lat.computeNtau();
    //lat.convergence();
    //lat.tau_o_plot();
    return 0;
     
}
