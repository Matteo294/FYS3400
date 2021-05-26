#include<iostream>
#include<fstream>
#include<cmath>

#define NSTEPS 1000
#define NSTEPS_NR 100

using namespace std;

double kB, Eg, Ed, Ec, ED;
double m_hole, m_elec, me;
int Nd;
double hbar;
double eV_to_Joule;

double ktExp(double x, double T){
    return exp(x/(kB*T));
}

double Ei(double T){
    return 0.5*Eg + (double)3/4*kB*T*log(m_hole/m_elec);
}

double ni(double T){
    return 2 * pow(kB*T/(2*M_PI*pow(hbar, 2)), 3/2) * pow(m_elec*m_hole, 3/4) * exp(-(Ec-ED)/(2*kB*T)) /1e6;
}

double n(double x, double T){
    return ni(T)*ktExp(x-Ei(T), T);
}

double p(double x, double T){
    return ni(T)*ktExp(Ei(T)-x, T);
}

double Ndplus(double x, double T, double ED){
    return Nd / (1 + 2*ktExp(x-ED, T));
}

/*double eq(double x, double T, double ED){
    return p(x, T) + Ndplus(x, T, ED) - n(x, T);
}

double eqdot(double x, double T, double ED){
    return - (double) 1/(kB*T) * (p(x, T) + n(x, T) + (double) 2/Nd*ktExp(x-ED, T)*pow(Ndplus(x, T, ED), 2));
}*/
double eq(double x, double T, double ED){
    return 2*ni(T)*sinh((x-Ei(T))/(kB*T)) + Ndplus(x, T, ED);
}
double eqdot(double x, double T, double ED){
    //cout << (kB*T) << " " << 1/pow(1+2*exp((x-ED)/(kB*T)), 2) << endl;
    return (double)2/(kB*T)*ni(T)*cosh((x-Ei(T))/(kB*T)) + (double)2/(kB*T)*Nd/pow(1+2*exp((x-ED)/(kB*T)), 2)*ktExp(x-Ed, T);
}



int main(){

    ofstream myfile;
    myfile.open("data.csv");
    myfile << "T,n" << endl;

    eV_to_Joule = 1.602e-19;
    kB = 1.381e-23/eV_to_Joule;
    Eg = 1.12; // eV
    Ed = 45e-3; // eV
    ED = Ec - Ed;
    Ec = Eg;
    me = 9.109e-31;
    m_elec = 1.06*me;
    m_hole = 0.56*me;
    Nd = (int) 1e17;
    hbar = 6.626e-34/(2*M_PI)/eV_to_Joule;
    
    double dT = (double) (1000-100)/NSTEPS;
    double dx;
    double T;

    for(T=100; T<1500; T+=dT){
        double x = 1;
        for(int i=0; i<NSTEPS_NR; i++){
            dx = (double) eq(x, T, ED)/eqdot(x, T, ED);
            cout << x << " " << x-Ei(T) << " " << eqdot(x, T, ED) << endl;
            x -= dx;
        }
        cout << endl;
        myfile << 1000/T << "," << n(x, T) << endl;
    }


    return 0;
}