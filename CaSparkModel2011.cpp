/* Monte-Carlo model of Ca spark from single cluster of RyRs
   MATLAB model written by Eric Sobie, Icahn School of Medicine, Mount Sinai
   Parameters specified in:
   Ramay et. al, Cardiovascular Research, 2011
   Written by: Ananya Pavuluri */

#include <iostream>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <time.h>

using namespace std; 

// FUNCTION DECLARATIONS
double gen_random();
bool is_a_timepoint(double time, double timearray[], int length);
double onedec(double x);

// ROUNDS A DOUBLE TO ONE DECIMAL PLACE
double onedec(double x)
{
      return round(x * 10) / 10;
}    

// CHECKS IF DATA AT A CERTAIN TIME SHOULD BE PLOTTED
// BY COMPARING TO PRE-SET ARRAY OF SIGNIFICANT TIMEPOINTS
bool is_a_timepoint(double time, double timearray[], int length) 
{
	for (int i = 0; i < length; i++) {
		
		if(onedec(timearray[i]) == time) {
			return true;
		} 
	}
	return false; 
}

// GENERATES A RANDOM DOUBLE BETWEEN 0 AND 1 FOR STOCHASTIC MODELING
double gen_random() 
{
	return double(rand()) / (double(RAND_MAX) + 1.0);
}

// VARIABLE DECLARATIONS AND SIMULATION
int main()
{

// random generator seed
srand (time(NULL));

// time steps
double dt = 0.00001000;
const double dt_record = 0.1000000;

// Open single RyR at interval, run for 'timeafter' ms
double interval = 0;
double timeafter = 30.0;

// number of sparks to simulate
const int trials = 3;

// Faraday's Constant. Used to convert flux to current
const double F = 96.485; // C/mmol

// Geometrical parameters
const double V_ss = 1.0000e-12;
const double V_JSR = 1.6000e-12; 

// Time constants, in ms
	// time constant for ...
const double tau_efflux = 1.78e-3;
	// time constant for NSR to JSR refilling
const double tau_refill = 6.5;

//Coupling energy between RyRs
const double EJequiv = 0.1;

// RyR permeability constant
const double D_ryr = 2.2e-12;

// RyR gating parameters
double kr_minus, kr_plus_max, Km_r_max, alpha_r, hill;
kr_minus = 0.48;        // max close rate, ms^-1
kr_plus_max = 30.0;         // max open rate, ms^-1
Km_r_max = 19.87;           // sensitivity of opening to subspace Ca, uM
alpha_r = 1.0e-3;           // luminal dependence factor 
hill = 4;                   // exponent

// # of RyR channels in a cluster
const int N_RyR = 28;

// Coupling rate
const double kcoup  = exp(2 * EJequiv / (N_RyR - 1));

// Subspace buffers 
	
	// Total [calmodulin], Total [SR membrane buffer],  Total [SL membrane buffer]
double bt[3] = {24, 47, 900}; // uM
	// On rates for calmodulin, SR membrane, SL membrane
double k_on[3] = {100, 115, 115};
	// Off rates for calmodulin, SR membrane, SL membrane
double k_off[3] = {38, 100, 1000};
	// converting from s^-1 to ms^-1
for(int n = 0; n < 3; n++) {
    k_on[n] *= 1e-3;
    k_off[n] *= 1e-3;
}

// JSR buffer calsequestrin
double CSQ = 30e3; // uM  // total [CSQ]
double KCSQ = 630; // uM  // Ca dissociation constant

// Fixed ionic concentrations
double Ca_myo = 0.1;        // bulk myoplasm [Ca2+]
double Ca_NSR = 1000;		  // NSR [Ca2+]
double Ca_ss = 0.1;         // subspace [Ca2+]
double Ca_JSR = 1000;       // JSR [Ca2+]

//Initializing arrays to hold results
double t_end = interval + timeafter;

int iterations = nearbyint(t_end/dt);
int outputs = nearbyint(t_end/dt_record) +1;

int plottime_size = (int) outputs * dt_record * 10;


// ALLOCATING HEAP MEMORY FOR DYNAMIC ARRAYS
	
	// times for the plot
double *plottime = new double [plottime_size];
double timepoint = 0;
for(int i = 0; i < plottime_size; i++) {
	plottime[i] = timepoint;
	timepoint = timepoint + dt_record;
}

	// Array for [Ca2+] in the JSR
double **Ca_JSR_all = new double * [outputs];
for (int i = 0; i < outputs; i++) {
		Ca_JSR_all[i] = new double [trials];
}

	// Array for subspace [Ca2+]
double **Ca_ss_all = new double * [outputs];
for (int i = 0; i < outputs; i++) {
		Ca_ss_all[i] = new double [trials];
}

	// Array for current released from RyR
double **Irel_all = new double * [outputs];
for (int i = 0; i < outputs; i++) {
		Irel_all[i] = new double [trials];
}

	// Array for number of receptors open
double **Nopen_all = new double * [outputs];
for (int i = 0; i < outputs; i++) {
		Nopen_all[i] = new double [trials];
}


// Initializing the arrays
for (int i = 0; i < outputs; i++) {
	for (int j = 0; j < trials; j++) {
		Ca_ss_all[i][j] = 0;
		Ca_JSR_all[i][j] = 0;
		Irel_all [i][j] = 0;
		Nopen_all[i][j] = 0;
	}
}

// Making array for buffers
double b[3];

double J_ryr = 0;
int count = 0;
				// ---------------------        SIMULATION       --------------------- // 
for (int i = 0; i < trials; i++) {
	
	Ca_ss = Ca_myo;
	Ca_JSR = Ca_NSR;
	double nopen = 0;
	
	// buffers
	for(int i = 0; i < 3; i++) {
	b[i] = (bt[i] * (k_off[i] / k_on[i]) ) /(k_off[i] / k_on[i] + Ca_ss) ;
	}

	int writedex = 0;
	double tlast = -1 * dt;
	int J_d = 0;
	bool neverspark = true;
	double time = 0.0;
	
	for (int j = 0; j < iterations; j++) {
		if (time >= interval && tlast < interval) {
			nopen = nopen + 5;
		}
		else if (time >= interval + 10 && nopen < 1 && neverspark)  {
			break;
		}
		double nclosed = N_RyR - nopen;
		
		// Fluxes and currents	
		J_ryr = nopen * D_ryr * (Ca_JSR - Ca_ss) / V_ss; // uM/ms
		double I_ryr = 1e6 * J_ryr * 2 * F * V_ss ; 	 // pA
		double J_efflux = (Ca_myo - Ca_ss) / tau_efflux ;
		double J_refill = (Ca_NSR - Ca_JSR) / tau_refill;

		// Buffers
		double db_dt[3];
		double J_buff = 0;
		for(int i = 0; i < 3; i++) {
			db_dt[i] = -1 * k_on[i] * b[i] * Ca_ss + k_off[i] * (bt[i] - b[i]);
			J_buff += db_dt[i];
		}

		double denom = pow(KCSQ + Ca_JSR, 2);
		double B_JSR = pow((1 + CSQ * KCSQ / denom), -1);

		// Writing arrays after the fluxes are calculated, before integration
		// and state switching

		if (is_a_timepoint(time, plottime, plottime_size)) {
				dt = 0.0000100; // the value of dt changes to 9.99999974737875e-06 
						   		// within this conditional for an unknown reason.
						   		// here, it is assigned 1e-5 again to prevent this.
				Irel_all[writedex][i] = I_ryr;
				Nopen_all[writedex][i] =  nopen;
				Ca_ss_all[writedex][i] = Ca_ss;
				Ca_JSR_all[writedex][i] = Ca_JSR;
				writedex = writedex + 1;

		}
		
		double Km_r = Km_r_max - alpha_r * Ca_JSR;
		double pow1 = pow(Ca_ss, hill);
		double pow2 =  pow(Km_r, hill);
		double kr_plus = kr_plus_max * pow1 / (pow1 + pow2);

		// Stochastic variables
		double  pincrease = dt * nclosed * kr_plus * pow(kcoup,2 * nopen + 1 - N_RyR);
		double  pdecrease = dt * nopen * kr_minus * pow(kcoup, 2 * nclosed + 1 - N_RyR);

		if (gen_random() < pincrease) {
			nopen += 1;
		}
		if (gen_random() < pdecrease) {
			nopen -= 1;
		}
		if(nopen >= 5) {
			neverspark = false;
		}

		// Subspace [Ca2+] and JSR [Ca2+] derivatives
		double dCass_dt = J_efflux + J_d + J_ryr + J_buff;
		double dCaJSR_dt = B_JSR * ( J_refill - J_ryr * V_ss/V_JSR);

		Ca_ss = Ca_ss + dt * dCass_dt;
		Ca_JSR = Ca_JSR + dt * dCaJSR_dt;

		// updating buffers
		for (int i = 0; i < 3; i++) {
			b[i] += dt*db_dt[i];
		}

		tlast = time;
		time += dt;
		// accounting for precision problems with floating point numbers.
		// ensures that time is accurate to 5 decimal places.
		time = round(time * 100000) / 100000;
	}

	// Write values at last time point t_end after the loop terminates
	Ca_ss_all [outputs-1][i]  = Ca_ss;
	Ca_JSR_all [outputs-1][i] = Ca_JSR;
	Irel_all[outputs-1][i]    = 1e6 * J_ryr * 2 * F * V_ss;
	Nopen_all[outputs-1][i]   = nopen;

}
	
	// WRITING OUTPUT TO CSV FILES

	ofstream Nopen;
	Nopen.open("N_open.csv");
	while (Nopen.is_open()) {
		for(int i = 0; i < outputs; i++) {
			Nopen << Nopen_all[i][0];
			for (int j = 1; j < trials; j++) {
				Nopen << "," << Nopen_all[i][j];
			}
			Nopen << endl;
		}
		break;
	}
	Nopen.close();
	ofstream Irel;
	Irel.open("Irel.csv");
	while (Irel.is_open()) {
		for(int i = 0; i < outputs; i++) {
			Irel << Irel_all[i][0];
			for (int j = 1; j < trials; j++) {
				Irel << "," << Irel_all[i][j];
			}
			Irel << endl;
		}
		break;
	}
	Irel.close();
	ofstream Cads;
	Cads.open("Ca_ss.csv");
	while (Cads.is_open()) {
		for(int i = 0; i < outputs; i++) {
			Cads << Ca_ss_all[i][0];
			for (int j = 1; j < trials; j++) {
				Cads << "," << Ca_ss_all[i][j];
			}
			Cads << endl;
		}
		break;
	}
	Cads.close();
	ofstream CaJSR;
	CaJSR.open("CaJSR.csv");
	while (CaJSR.is_open()) {
		for(int i = 0; i < outputs; i++) {
			CaJSR << Ca_JSR_all[i][0];
			for (int j = 1; j < trials; j++) {
				CaJSR << "," << Ca_JSR_all[i][j];
			}
			CaJSR << endl;
		}
		break;
	}
	ofstream plot_time;
	plot_time.open("plottime.csv");
	while (plot_time.is_open()) {
		for(int i = 0; i < outputs; i++) {
			plot_time << plottime[i];
			if (i != outputs-1) 
				plot_time << ",";
		}
		break;
	}
	plot_time.close(); 
	// Recycling memory
	for (int k = 0; k < trials; k++) {
		delete [] Ca_ss_all[k];
		delete [] Ca_JSR_all[k];
		delete [] Nopen_all[k];
		delete [] Irel_all[k];
	}
	delete [] plottime;
	delete [] Ca_ss_all;
	delete [] Ca_JSR_all;
	delete [] Nopen_all;
	delete [] Irel_all;
	
	cout << "End of program." << endl;

	return 0;
}


