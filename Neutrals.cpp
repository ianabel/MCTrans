
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

double IonNeutralCrossSection( double Ti )
{
	// Need to implement cross section dependent on the FuelName
	return 0.6937E-14 * ::pow( 1 - 0.155 * ::log10( Ti ), 2 ) / (1 + 0.1112E-14 * ::pow( Ti, 3.3 ));
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
