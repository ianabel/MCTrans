
#include "MirrorPlasma.hpp"

#include "PlasmaPhysics.hpp"

#include <cmath>
#include <iostream>
#include <vector>
#include <string>

#include <boost/math/quadrature/trapezoidal.hpp>

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
	return density * sigmaV;
}

// Fits for sigma are given starting on pg 234 in Janev, 1987
// <cite>
// Energy in electron Volts

constexpr unsigned int N_JANEV_COEFFS = 9;

double EvaluateJanevCrossSectionFit ( std::vector<double> PolynomialCoefficients, double Energy )
{
	if ( PolynomialCoefficients.size() != N_JANEV_COEFFS )
		throw std::invalid_argument( "Janev uses fixed order fits, there are not " + std::to_string( N_JANEV_COEFFS ) + " numbers. Something is wrong." );
	double sum = 0.0;
	for ( size_t n = 0; n < PolynomialCoefficients.size(); n++ ) {
      sum += PolynomialCoefficients.at( n ) * ::pow( ::log(Energy), n );
   }
	double sigma = ::exp(sum);
	return sigma;
}

class CrossSection {
	public:
		CrossSection( std::function<double( double )> sigma, double E_min, double E_max, double mu, double mr )
		{
			sigmaImplementation = sigma;
			MinEnergy = E_min;
			MaxEnergy = E_max;
			ReducedMass = mu;
			RelativeMass = mr;
		};
		
		// Masses in units of the proton mass
		double ReducedMass,RelativeMass;
		// in units of electron volts (not keV!)
		double MinEnergy,MaxEnergy;

		// Return the cross section in cm^2 -- not barns, not m^2
		double operator()( double CentreOfMassEnergy ) const
		{
			return sigmaImplementation( CentreOfMassEnergy );
		}

	private:
		std::function<double( double )> sigmaImplementation;
};


double IonNeutralCrossSection( double Ti )
{
	// Need to implement cross section dependent on the FuelName
	// Janev 1987 3.1.8
	double sigma_1s = 0.6937e-14 * ::pow( 1 - 0.155 * ::log10( Ti ), 2 ) / (1 + 0.1112E-14 * ::pow( Ti, 3.3 ));

	// Janev 1987 3.1.9
	std::vector<double> aSigma_2p = {-2.197571949935e+01, -4.742502251260e+01, 3.628013140596e+01, -1.423003075866e+01, 3.273090240144e+00, -4.557928912260e-01, 3.773588347458e-02, -1.707904867106e-03, 3.251203344615e-05};
	// Janev 1987 3.1.10
	std::vector<double> aSigma_2s = {-1.327325087764e+04, 1.317576614520e+04, -5.683932157858e+03, 1.386309780149e+03, -2.089794561307e+02, 1.992976245274e+01, -1.173800576157e+00, 3.902422810767e-02, -5.606240339932e-04};

	// Contribution from ground -> 2p orbital
	double sigma_2p;
	if ( Ti < 19.0 ) {
		sigma_2p = 0;
	} else {
		sigma_2p = EvaluateJanevCrossSectionFit( aSigma_2p, Ti );
	}

	// Contribution from ground -> 2s orbital
	double sigma_2s;
	if ( Ti < 0.1 ) {
		sigma_2s = 0;
	} else {
		sigma_2s = EvaluateJanevCrossSectionFit( aSigma_2s, Ti );
	}

	return sigma_1s + sigma_2p + sigma_2s;
}



double rateCoeff( double Ti, CrossSection const & sigma )
{
	// E and T in eV, sigma in cm^2
	// k=<σv> in m^3/s
	// mu, mr is in proton masses
	// Assumes a Maxwellian distribution, COM energy (as opposed to incident)
	Ti *= ElectronCharge; // Convert to Joules

	auto integrand = [&]( double Energy ) {
		double sigmaM2 = sigma( Energy ) * 1e-4; // sigma is in cm^2, we need m^2
		return Energy * sigmaM2 * ::exp( -Energy / Ti );
	};

	constexpr double tolerance = 1e-6;
	return ( 4.0 / ( ::sqrt(2 * M_PI * sigma.ReducedMass * ProtonMass * Ti) * Ti ) ) * boost::math::quadrature::trapezoidal( integrand, sigma.MinEnergy, sigma.MaxEnergy, tolerance );
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
