#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <hdf5.h>
#include <mpi.h>
#include <stdbool.h>
#include <unistd.h>
#include <math.h>


double* calcHydrogenNumberDensity(double** gas_metallicities, double* gas_densities, int Ngas)
{
    double *nH; 
    nH = (double *) malloc (Ngas * sizeof (double));
    int i;
    double h_massfrac;
    double unit_M = pow(10,10) * 1.98855 * pow(10,33);
    double proton_mass = 1.6726219*pow(10,-27)*(1000.0/unit_M); //appropriate units
    for (i=0;i<Ngas;i++) {
        h_massfrac = 1 - (gas_metallicities[i][0]+gas_metallicities[i][1]);
        nH[i] = (gas_densities[i]*h_massfrac) / proton_mass;
    }
    return nH;
}


double* calcH1Abundance(double* gas_masses, double* neutral_hydrogen_densities, double* kernalLengths, double* gas_densities, double** gas_metallicities, int Ngas)
{
    int i;
    double Z_MW = 0.02; //Assuming Milky way metallicity ~ solar on average
    double Z_solar = 0.02; //From Gizmo Documentation
    double unit_M = pow(10,10)*1.98855*pow(10,33);
    double unit_L = 3.086*pow(10,21);
    double mu_H = 2.3*pow(10,-27)*1000/unit_M;
    double M_H = 1.67353269159582103*pow(10,-27)*1000/unit_M;

    double Z,sobColDens,tau,chi,s,fH1;
    double *fH2, *NH1;
    fH2 = (double *) malloc (Ngas * sizeof (double));
    NH1 = (double *) malloc (Ngas * sizeof (double));
    double epsilon = pow(10,-20);
    for (i=0;i<Ngas;i++) {
        Z = gas_metallicities[i][0];
        sobColDens=kernalLengths[i]*gas_densities[i];
        tau = sobColDens*Z / (mu_H*Z_MW) * pow(10,-21) / pow(unit_L,2);
        if (tau==0){tau=epsilon;}//Avoid Divide by zero
        
        chi = 3.1 * (1+3.1*pow(Z/Z_solar,0.365)) / 4.1; //Approximation

        s = log(1+0.6*chi+0.01*pow(chi,2)) / (0.6*tau);
        if (s==-4){s==-4+epsilon;}//Avoid Divide by zero
        
        fH2[i] = (1-0.5*s) / (1+0.25*s);
        if (fH2<0){fH2=0;}//Get rid of nonphysical negative values
        fH1 = neutral_hydrogen_densities[i] * (1-Z-gas_metallicities[i][1]) * (1-fH2[i]);
        NH1[i] = gas_masses[i]*fH1 / M_H;

        
    }

    return NH1;

}

double* calcTemperatures(double* internalEnergy,double* electronAbundance,double** metallicity, int Npart)
{
  int i;
  double unit_L = 3.086*pow(10,21);
  double unit_V = 1.0*pow(10,5);
  double unit_T = unit_L/unit_V;
  double unit_M = pow(10,10) * 1.98855 * pow(10,33);
  double gamma = 5/3;
  double k_Boltzmann = 1.38064852*pow(10,-23)*pow((100.0/unit_L),2)*(1000.0/unit_M)*pow(unit_T,2);
  double proton_mass = 1.6726219*pow(10,-27)*(1000.0/unit_M);
  double mu,mean_molecular_weight,helium_mass_fraction,y_helium;
  double* temperatures;
  temperatures = (double *) malloc (Npart * sizeof (double));
  for (i=0;i<Npart;i++) {
    helium_mass_fraction = metallicity[i][1];
    y_helium = helium_mass_fraction / (4*(1-helium_mass_fraction));
    mu = (1 + 4*y_helium) / (1+y_helium+electronAbundance[i]);
    mean_molecular_weight = mu*proton_mass;
    temperatures[i] = mean_molecular_weight * (gamma-1) * internalEnergy[i] / k_Boltzmann;
  }

return temperatures;
}

