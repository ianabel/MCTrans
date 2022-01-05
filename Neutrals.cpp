
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

double hydrogenIonizationFraction( double ElectronDensity, double ElectronTemperature)
{
	// Analytical solution from the equilibrium solution to Saha equation
	// derivation here: http://www.astro.wisc.edu/~townsend/resource/teaching/astro-310-F09/hydrogen-ionization.pdf
	constexpr double hydrogenIonizationPotential = 13.6; // in units of eV
	return ::pow( 2.0 * M_PI * ElectronMass * ElectronTemperature / ::pow( PlanckConstant, 2.0 ), 1.5 ) * ::exp( -hydrogenIonizationPotential / ElectronTemperature ) / ElectronDensity;
}

bool isShortMeanFreePathRegime( double meanFreePath, double characteristicLength, double minRatio = 10.0 )
{
	return ( characteristicLength / meanFreePath ) >= minRatio;
}

bool isCoronalEquilibrium( double excitationRate, double deexcitationRate, double minRatio = 10.0 )
{
	return ( deexcitationRate / excitationRate ) >= minRatio;
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

	// Minimum energy of cross section in eV
	const double minimumEnergySigma_1s = 0.1;
	const double minimumEnergySigma_2p = 19.0;
	const double minimumEnergySigma_2s = 0.1;

	// Janev 1987 3.1.9
	// p + H(1s) --> H(2p) + p
	std::vector<double> aSigma_2p = {-2.197571949935e+01, -4.742502251260e+01, 3.628013140596e+01, -1.423003075866e+01, 3.273090240144e+00, -4.557928912260e-01, 3.773588347458e-02, -1.707904867106e-03, 3.251203344615e-05};
	// Janev 1987 3.1.10
	// p + H(1s) --> H(2s) + p
	std::vector<double> aSigma_2s = {-1.327325087764e+04, 1.317576614520e+04, -5.683932157858e+03, 1.386309780149e+03, -2.089794561307e+02, 1.992976245274e+01, -1.173800576157e+00, 3.902422810767e-02, -5.606240339932e-04};

	// Contribution from ground -> ground state
	// Janev 1987 3.1.8
	// p + H(1s) --> H(1s) + p
	double sigma_1s;
	if ( Ti < minimumEnergySigma_2p ) {
		sigma_1s = 0;
	} else {
		sigma_1s = 0.6937e-14 * ::pow( 1 - 0.155 * ::log10( Ti ), 2 ) / (1 + 0.1112e-14 * ::pow( Ti, 3.3 ));
	}

	// Contribution from ground -> 2p orbital
	double sigma_2p;
	if ( Ti < minimumEnergySigma_2p ) {
		sigma_2p = 0;
	} else {
		sigma_2p = EvaluateJanevCrossSectionFit( aSigma_2p, Ti );
	}

	// Contribution from ground -> 2s orbital
	double sigma_2s;
	if ( Ti < minimumEnergySigma_2s ) {
		sigma_2s = 0;
	} else {
		sigma_2s = EvaluateJanevCrossSectionFit( aSigma_2s, Ti );
	}

	return sigma_1s + sigma_2p + sigma_2s;
}

double electronHydrogenIonizationCrossSection( double Te )
{
	// Minimum energy of cross section in eV
	const double ionizationEnergy = 13.6;
	const double minimumEnergySigma = ionizationEnergy;

	// Contribution from ground state
	// Janev 1993, ATOMIC AND PLASMA-MATERIAL INTERACTION DATA FOR FUSION, Volume 4
	// Equation 1.2.1
	// e + H(1s) --> e + H+ + e
	// Accuracy is 10% or better
	double fittingParamA = 0.18450;
	std::vector<double> fittingParamB = { -0.032226, -0.034539, 1.4003, -2.8115, 2.2986 };

	double sigma;
	if ( Te < minimumEnergySigma ) {
		sigma = 0;
	}
	else {
		double sum = 0.0;
		for ( size_t n = 0; n < fittingParamB.size(); n++ ) {
	      sum += fittingParamB.at( n ) * ::pow( 1 - ionizationEnergy / Te, n );
	   }
		sigma = 1.0e-13 / ( ionizationEnergy * Te ) * ( fittingParamA * ::log( Te / ionizationEnergy ) + sum );
	}
	return sigma;
}

double evaluateJanevDFunction( double beta )
{
	// Function from Janev 1987 Appendix C
	// Used in some analytical fits for cross sections
	double DFunctionVal;
	if ( beta < 1e-3 ){
		DFunctionVal = 4 * beta * ::log( 1.4 / beta );
	} else if ( beta > 10 ) {
		DFunctionVal = beta / 2 * ::exp( - ::sqrt( 2 * beta ) );
	} else {
		// Create look up table and linearly interpolate to find DFunction
		std::vector<double> betaInverse = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 8.0, 10.0, 12.5, 15.0, 20.0, 25.0, 30.0, 40.0, 60.0, 80.0, 100.0, 15.0, 200.0, 300.0, 400.0, 600.0, 800.0, 1000.0 };
		std::vector<double> DFunction = { 0.057, 0.104, 0.135, 0.157, 0.175, 0.194, 0.230, 0.264, 0.296, 0.328, 0.367, 0.388, 0.399, 0.405, 0.410, 0.405, 0.399, 0.380, 0.345, 0.318, 0.285, 0.263, 0.227, 0.205, 0.190, 0.168, 0.141, 0.124, 0.110, 0.092, 0.080, 0.054, 0.042, 0.035, 0.028 };
		for (size_t i = 0; i < betaInverse.size(); i++ ){
			if ( 1.0 / beta < betaInverse.at( i ) ){}
			else {
				double x1 = betaInverse.at( i - 1 );
				double x2 = betaInverse.at ( i );
				double y1 = DFunction.at( i - 1 );
				double y2 = DFunction.at ( i );
				DFunctionVal = y1 + ( y2 - y1 ) * ( 1.0 / beta - x1 ) / ( x2 - x1 );
			}
		}
	}
	return DFunctionVal;
}

double protonHydrogenIonizationCrossSection1( double Ti )
{
	// Minimum energy of cross section in eV
	const double minimumEnergySigma = 13.6;

	// Contribution from ground state
	// Janev 1987 3.1.6
	// p + H(1s) --> p + H+ + e
	const double lambda_eff = 0.808;
	const double lambda_01 = 0.7448;
	const double omega_i = 0.5;
	const double omega_01 = 0.375;
	double v = 6.3246e-3 * ::sqrt( Ti );
	double beta_i = lambda_eff * omega_i / ::pow( v, 2 );
	double beta_01 = lambda_01 * omega_01 / ::pow( v, 2 );

	double sigma;
	if ( Ti < minimumEnergySigma ) {
		sigma = 0;
	} else {
		sigma = 1.76e-16 * ( lambda_eff * evaluateJanevDFunction( beta_i ) / omega_i  + lambda_01 * evaluateJanevDFunction( beta_01 ) / ( 8 * omega_01 ) );
	}
	return sigma;
}

double protonHydrogenIonizationCrossSection2( double Ti )
{
	// Minimum energy of cross section in keV
	const double minimumEnergySigma = 0.2;
	double TiKEV = Ti / 1000;

	// Contribution from ground state
	// Janev 1993, ATOMIC AND PLASMA-MATERIAL INTERACTION DATA FOR FUSION, Volume 4
	// Equation 2.2.1
	// H+ + H(1s) --> H+ + H+ + e
	// Accuracy is 30% or better
	const double A1 = 12.899;
	const double A2 = 61.897;
	const double A3 = 9.2731e3;
	const double A4 = 4.9749e-4;
	const double A5 = 3.9890e-2;
	const double A6 = -1.5900;
	const double A7 = 3.1834;
	const double A8 = -3.7154;

	double sigma;
	if ( TiKEV < minimumEnergySigma ) {
		sigma = 0;
	}
	else {
		// Energy is in units of keV
		sigma = 1e-16 * A1 * ( ::exp( -A2 / TiKEV ) * ::log( 1 + A3 * TiKEV ) / TiKEV + A4 * ::exp( -A5 * TiKEV ) / ( ::pow( TiKEV, A6 ) + A7 * ::pow( TiKEV, A8 ) ) );
	}
	return sigma;
}

double electronHydrogenExcitationN2CrossSection( double Te )
{
	// Minimum energy of cross section in eV
	const double excitationEnergy = 10.2;
	const double cutoffEnergy1 = 11.56;
	const double cutoffEnergy2 = 12.23;

	// Contribution from ground state
	// Janev 1993, ATOMIC AND PLASMA-MATERIAL INTERACTION DATA FOR FUSION, Volume 4
	// Equation 1.1.3
	// e + H(1s) --> e + H*(n=2)
	// Accuracy is 10% or better
	std::vector<double> fittingParamA = { 1.4182, -20.877, 49.735, -46.249, 17.442, 4.4979 };

	double sigma;
	if ( Te < excitationEnergy ) {
		sigma = 0;
	}
	else if ( Te < cutoffEnergy1 ) {
		sigma = 1e-16 * ( 0.255 + 0.1865 * ( Te - excitationEnergy ) );
	}
	else if ( Te < cutoffEnergy2 ) {
		sigma = 5.025e-17;
	}
	else {
		double sum = 0.0;
		double XEnergy = Te / excitationEnergy;
		for ( size_t n = 0; n < fittingParamA.size() - 1; n++ ) {
	      sum += fittingParamA.at( n ) / ::pow( XEnergy, n - 1 );
	   }
		sigma = 5.984e-16 / Te * ( sum + fittingParamA.back() * ::log( XEnergy ) );
	}
	return sigma;
}

double protonHydrogenExcitationN2CrossSection( double Ti )
{
	// Minimum energy of cross section in keV
	const double minimumEnergySigma = 0.6;
	double TiKEV = Ti / 1000;

	// Contribution from ground state
	// Janev 1993, ATOMIC AND PLASMA-MATERIAL INTERACTION DATA FOR FUSION, Volume 4
	// Equation 2.2.1
	// H+ + H(1s) --> H+ + H+ + e
	// Accuracy is 100% or better
	const double A1 = 34.433;
	const double A2 = 44.057;
	const double A3 = 0.56870;
	const double A4 = 8.5476;
	const double A5 = 7.8501;
	const double A6 = -9.2217;
	const double A7 = 1.8020e-2;
	const double A8 = 1.6931;
	const double A9 = 1.9422e-3;
	const double A10 = 2.9068;

	double sigma;
	if ( TiKEV < minimumEnergySigma ) {
		sigma = 0;
	}
	else {
		// Energy is in units of keV
		sigma = 1e-16 * A1 * ( ::exp( -A2 / TiKEV ) * ::log( 1 + A3 * TiKEV ) / TiKEV + A4 * ::exp( -A5 * TiKEV ) / ( ::pow( TiKEV, A6 ) ) + A7 * ::exp( -A8 / TiKEV ) / ( 1 + A9 * ::pow( TiKEV, A10 ) ) );
	}
	return sigma;
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


// MFP of hydrogen atom in the plasma
double NeutralCXMeanFreePath( double IonDensity, double IonTemperature )
{
	// TODO
	throw std::logic_error( "Unimplemented Function." );
	return 0.0;
}

// MFP of proton in neutral atomic hydrogen
double IonCXMeanFreePath( double NeutralDensity, double IonTemperature )
{
	// TODO
	throw std::logic_error( "Unimplemented Function." );
	return 0.0;
}
