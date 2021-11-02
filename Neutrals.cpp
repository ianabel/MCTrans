
#include "MirrorPlasma.hpp"
#include <cmath>
#include <iostream>

constexpr double KB = 8.617E-5 // [eV/K]
constexpr double MP = 1.673E-27 // [kg]
constexpr double QE = 1.602e-19 // [C]

// Reference Cross-section for Charge-Exchange neutral cross-section
// in cm^2
constexpr double BaseCXCrossSection = 1e-15;
constexpr double BaseIonizationCrossSection = 1e-10;

double meanFreePath( double density, double crossSection )
{
	return 1.0 / ( density * crossSection );
}

// Density should be neutral density
double CXRate( double density, double sigmaV )
{
	return density * sigmaV
}

double trapz( double yi[], double xi[] )
{
	if ( xi.size() != yi.size() )
		throw std::invalid_argument( "The xi and yi lists must be the same length" );

	double sum = 0.0;
	for(int i = 0; i<xi.size(); i++){
      sum += 0.5 * ( yi[i] + yi[i + 1] ) / ( xi[i] - xi[i + 1] );
   }
	 return sum;
}

// Fits for sigma are given starting on pg 234 in Janev, 1987
double JanevCrossSectionFit ( double a[], double Energy )
{
	double sum = 0.0;
	for(int n = 0; n < an.size(); n++){
      sum += a[n] * ::pow( ::log(Energy), n )
   }
	double sigma = ::exp(sum)
	return sigma
}

double IonNeutralCrossSection( double Ti )
{
	// Need to implement cross section dependent on the FuelName
	// Janev 1987 3.1.8
	double sigma_1s = 0.6937E-14 * ::pow( 1 - 0.155 * ::log10( Ti ), 2 ) / (1 + 0.1112E-14 * ::pow( Ti, 3.3 ));

	// Janev 1987 3.1.9
	double aSigma_2p[9] = {-2.197571949935e+01, -4.742502251260e+01, 3.628013140596e+01,
												 -1.423003075866e+01, 3.273090240144e+00, -4.557928912260e-01,
												 3.773588347458e-02, -1.707904867106e-03, 3.251203344615e-05};
	if Ti < 19.0{
		double sigma_2p = 0
	} else{
		double sigma_2p = JanevCrossSectionFit(aSigma_2s, Ti)
	}

	// Janev 1987 3.1.10
	double aSigma_2s[9] = {-1.327325087764e+04, 1.317576614520e+04, -5.683932157858e+03,
												 1.386309780149e+03, -2.089794561307e+02, 1.992976245274e+01,
											 	 -1.173800576157e+00, 3.902422810767e-02, -5.606240339932e-04};
	if Ti < 0.1{
		double sigma_2s = 0
	} else{
		double sigma_2s = JanevCrossSectionFit(aSigma_2s, Ti)
	}

	return sigma_1s + sigma_2p + sigma_2s;
}

double rateCoeff( double Energy, double Ti, double sigma, double mu, double mr )
{
	// E and T in eV, sigma in cm^2
	// k=<σv> in m^3/s
	// mu, mr is in proton masses
	// Assumes a Maxwellian distribution, COM energy (as opposed to incident)
	sigma *= 1e-4; // Convert to m^2
	Energy *= QE * mr; // Convert to J and COM energy
	Ti *= QE; // Convert to J
	mu *= MP;
	return 4 / np.sqrt(2 * M_PI * mu * Ti) / Ti * trapz( Energy * sigma * ::exp( -Energy / Ti ), Energy );
}

// Ionization Rate as a function of
// electron density and electron temperature
// returned in ADAS-like units of m^3/s (so multiply by
// n_neutral * n_e * Volume to get rate)
double IonizationRate( double Ne, double Te )
{
	return 1.0;
}

// CX Rate as a function of
// ion density and ion temperature
double CXCrossSection( double Ni, double Ti )
{
	double IonThermalSpeed = ::sqrt( 2.0 * Ti * ReferenceTemperature / ElectronMass );
	return BaseCXCrossSection * 1e-4 * IonThermalSpeed;
}

double RecombinationRate( double Ne, double Te )
{
	return 0.0;
}

void MirrorPlasma::ComputeSteadyStateNeutrals()
{
	double CurrentIonizationRate = IonizationRate( ElectronDensity * ReferenceDensity, ElectronTemperature * ReferenceTemperature );
	// Assume mix of neutrals is such that QN is maintained so we just need to produce enough electrons from the source gas
	double ElectronLossRate = ParallelElectronParticleLoss() + ClassicalElectronParticleLosses();
	// Steady State requires Losses = n_N * n_e * IR
	// so n_N = Losses / ( n_e * IR )
	NeutralDensity = ElectronLossRate / ( ElectronDensity * ReferenceDensity * CurrentIonizationRate );
	// Normalize neutral density
	NeutralDensity = NeutralDensity / ReferenceDensity;
	NeutralSource = ElectronLossRate;
}

/*
 * To the neutrals the plasma is an onrushing torrent. Moving to the plasma rest fram, the neutrals are a high-energy beam ( width is
 * determined by the temperature of the neutral gas, which is low ).
 *
 * We can thus use beam stopping coefficients to try and work out the ionization rates..
 */
