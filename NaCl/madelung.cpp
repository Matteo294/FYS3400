#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

int main(){

    ofstream outfile("madelung.csv");

    // Energy conversion factor
    const double eV_to_Joule = 1.602e-19;
    // Factors to convert to atomic units
    const double E0 = 4.360e-18;
    const double x0 = 5.292e-11;
    // Lattice parameters
    const int N = 50; // Number of molecules
    const int a = 2.82e-10/x0; // Lattice constant
    const double lam = 1e3*eV_to_Joule/E0;
    const double rho = 0.32e-10;

    outfile << "N_atoms,Madelung,Madelung_borders" << endl;

    double s, s2; // Sums
    for(int n=2; n<=2*N; n+=2){

        // Method 1: takes account of the border effects
        s = 0.0;
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
                if(i != j){
                    s += (double) pow(-1, i+j)/(abs(i-j)*a);
                }
            }
        }

        // Method 2: consider symmetrical atoms
        s2 = 0.0;
        for(int j=1; j<(int)n/2; j++){
            s2 += (double) pow(-1, j)/(j*a);
        }

        // Write to the file (plots are made with python)
        outfile << n << "," << -2*s2*a << "," << -s*a/n << endl;
    }

    return 0;
}