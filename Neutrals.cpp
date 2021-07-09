
#include "MirrorPlasma.hpp"
#include <cmath>
#include <iostream>

// Reference Cross-section for Charge-Exchange neutral cross-section
// in cm^2
constexpr double BaseCXCrossSection = 1e-15;
constexpr double BaseIonizationCrossSection = 1e-10;

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
