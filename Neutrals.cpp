
#include "MirrorPlasma.hpp"
#include "PlasmaPhysics.hpp"
#include "AtomicPhysics.hpp"

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <boost/math/quadrature/gauss_kronrod.hpp>

double neutralsRateCoefficientHot( CrossSection const & sigma, MirrorPlasma const & plasma )
{
	// E and T in eV, sigma in cm^2
	// k=<σv> in m^3/s
	// Assumes a Maxwellian distribution, COM energy (as opposed to incident)

	double temperature;
	if ( sigma.Particle.Name == "Electron" ){
		temperature = plasma.ElectronTemperature * ReferenceTemperature; // Convert to Joules
	}
	else{
		temperature = plasma.IonTemperature * ReferenceTemperature; // Convert to Joules
	}

	double Jacobian = ElectronCharge; // The integral is over Energy, which is in units of electronvolts, so transform the integrand back to eV
	auto integrand = [&]( double Energy ) {
		double sigmaM2 = sigma( Energy ) * 1e-4; // sigma is in cm^2, we need m^2
		return Energy * ElectronCharge * sigmaM2 * ::exp( -Energy * ElectronCharge / temperature ) * Jacobian;
	};

	constexpr double tolerance = 1e-5;
	constexpr unsigned MaxDepth = 10;
	double HotRateCoeff = 4.0 / ( ::sqrt(2 * M_PI * sigma.ReducedMass * temperature) * temperature )
	         * boost::math::quadrature::gauss_kronrod<double, 255>::integrate( integrand, sigma.MinEnergy, sigma.MaxEnergy, MaxDepth, tolerance );

#if defined( DEBUG ) && defined( ATOMIC_PHYSICS_DEBUG )
	std::cerr << "Computing a cold rate coefficient at T = " << plasma.ElectronTemperature/1000 << " eV and M = " << plasma.MachNumber << " gave <sigma v> = " << ColdRateCoeff  << std::endl;
#endif

	return HotRateCoeff;

}

double neutralsRateCoefficientCold( CrossSection const & sigma, MirrorPlasma const & plasma )
{
	// sigma in cm^2
	// k=<σv> in m^3/s
	// Assumes ions are Maxwellian and neutrals are stationary
	double temperature;
	if ( sigma.Particle.Name == "Electron" ){
		temperature = plasma.ElectronTemperature * ReferenceTemperature; // Convert to Joules
	}
	else{
		temperature = plasma.IonTemperature * ReferenceTemperature; // Convert to Joules
	}

	double thermalSpeed = ::sqrt( 2.0 * temperature / sigma.Particle.Mass );
	double thermalMachNumber = plasma.MachNumber * ::sqrt( ::abs( sigma.Particle.Charge ) * plasma.ElectronTemperature * ReferenceTemperature / ( 2 * temperature ) );

	double Jacobian = ElectronCharge / (sigma.ReducedMass * thermalSpeed * thermalSpeed); // The integral is over Energy, which is in units of electronvolts, so transform the integrand back to eV, including change of variables from du to dE (less one power of u, which cancels with one in the integrand
	auto integrand = [&]( double Energy ) {
		double velocity = ::sqrt( 2.0 * Energy * ElectronCharge / sigma.ReducedMass );
		double u = velocity / thermalSpeed;
		double sigmaM2 = sigma( Energy ) * 1e-4; // sigma is in cm^2, we need m^2
		return u * sigmaM2 * ( ::exp( -::pow( thermalMachNumber - u, 2 ) ) - ::exp( -::pow( thermalMachNumber + u, 2 ) ) ) * Jacobian;
	};

	constexpr double tolerance = 1e-5;
	constexpr unsigned MaxDepth = 10;
	double ColdRateCoeff = thermalSpeed / ( thermalMachNumber * ::sqrt(M_PI) )
	        * boost::math::quadrature::gauss_kronrod<double, 255>::integrate( integrand, sigma.MinEnergy, sigma.MaxEnergy, MaxDepth, tolerance );

#if defined( DEBUG ) && defined( ATOMIC_PHYSICS_DEBUG )
	std::cerr << "Computing a cold rate coefficient at T = " << plasma.ElectronTemperature/1000 << " eV and M = " << plasma.MachNumber << " gave <sigma v> = " << ColdRateCoeff  << std::endl;
#endif

	return ColdRateCoeff;
}

/*
double meanFreePath( double density, double crossSection )
{
	return 1.0 / ( density * crossSection );
}

bool isShortMeanFreePathRegime( double meanFreePath, double characteristicLength, double minRatio = 10.0 )
{
	return ( characteristicLength / meanFreePath ) >= minRatio;
}

bool isCoronalEquilibrium( double excitationRate, double deexcitationRate, double minRatio = 10.0 )
{
	return ( deexcitationRate / excitationRate ) >= minRatio;
}
*/

double reactionRate( double density_1, double density_2, double rateCoefficient, std::shared_ptr<MirrorPlasma> pMirrorPlasma )
{
	return density_1 * density_2 * rateCoefficient * pMirrorPlasma->pVacuumConfig->PlasmaVolume();
}

double EvaluateJanevCrossSectionFit ( std::vector<double> PolynomialCoefficients, double Energy )
{
	constexpr unsigned int N_JANEV_COEFFS = 9;
	if ( PolynomialCoefficients.size() != N_JANEV_COEFFS )
		throw std::invalid_argument( "Janev uses fixed order fits, there are not " + std::to_string( N_JANEV_COEFFS ) + " numbers. Something is wrong." );
	double sum = 0.0;
	for ( size_t n = 0; n < PolynomialCoefficients.size(); n++ ) {
      sum += PolynomialCoefficients.at( n ) * ::pow( ::log(Energy), n );
   }
	double sigma = ::exp(sum);
	return sigma;
}

/*
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
*/

double electronImpactIonizationCrossSection( double CoMEnergy )
{
	// Minimum energy of cross section in eV
	constexpr double ionizationEnergy = 13.6;
	constexpr double minimumEnergySigma = ionizationEnergy;

	// Contribution from ground state
	// Janev 1993, ATOMIC AND PLASMA-MATERIAL INTERACTION DATA FOR FUSION, Volume 4
	// Equation 1.2.1
	// e + H(1s) --> e + H+ + e
	// Accuracy is 10% or better
	constexpr double fittingParamA = 0.18450;
	constexpr std::array<double,5> fittingParamB{ -0.032226, -0.034539, 1.4003, -2.8115, 2.2986 };

	double sigma;
	if ( CoMEnergy < minimumEnergySigma ) {
		sigma = 0;
	}
	else {
		double sum = 0.0;
		double x = 1.0 - ionizationEnergy / CoMEnergy;
		for ( size_t n = 0; n < fittingParamB.size(); n++ ) {
	      sum += fittingParamB.at( n ) * ::pow( x, n );
	   }
		sigma = 1.0e-13 / ( ionizationEnergy * CoMEnergy ) * ( fittingParamA * ::log( CoMEnergy / ionizationEnergy ) + sum );
	}
	return sigma;
}

// Energy in electron volts, returns cross section in cm^2
double protonImpactIonizationCrossSection( double Energy )
{
	// Minimum energy of cross section in keV
	const double minimumEnergySigma = 0.2;
	// Convert to keV
	double CoMEnergy = Energy / 1000;

	// Contribution from ground state
	// Janev 1993, ATOMIC AND PLASMA-MATERIAL INTERACTION DATA FOR FUSION, Volume 4
	// Equation 2.2.1
	// H+ + H(1s) --> H+ + H+ + e
	// Accuracy is 30% or better
	constexpr double A1 = 12.899;
	constexpr double A2 = 61.897;
	constexpr double A3 = 9.2731e3;
	constexpr double A4 = 4.9749e-4;
	constexpr double A5 = 3.9890e-2;
	constexpr double A6 = -1.5900;
	constexpr double A7 = 3.1834;
	constexpr double A8 = -3.7154;

	double sigma;
	if ( CoMEnergy < minimumEnergySigma ) {
		sigma = 0;
	}
	else {
		// Energy is in units of keV
		sigma = 1e-16 * A1 * ( ::exp( -A2 / CoMEnergy ) * ::log( 1 + A3 * CoMEnergy ) / CoMEnergy + A4 * ::exp( -A5 * CoMEnergy ) / ( ::pow( CoMEnergy, A6 ) + A7 * ::pow( CoMEnergy, A8 ) ) );
	}
	return sigma;
}

double HydrogenChargeExchangeCrossSection( double CoMEnergy )
{
	// Minimum energy of cross section in eV
	const double minimumEnergySigma_1s = 0.1;
	const double minimumEnergySigma_2p = 19.0;
	const double minimumEnergySigma_2s = 0.1;

	// Contribution from ground -> ground state
	// Janev 1987 3.1.8
	// p + H(1s) --> H(1s) + p
	double sigma_1s;
	if ( CoMEnergy < minimumEnergySigma_2p ) {
		sigma_1s = 0;
	} else {
		sigma_1s = 0.6937e-14 * ::pow( 1 - 0.155 * ::log10( CoMEnergy ), 2 ) / (1 + 0.1112e-14 * ::pow( CoMEnergy, 3.3 ));
	}

	// Janev 1987 3.1.9
	// p + H(1s) --> H(2p) + p
	std::vector<double> aSigma_2p = {-2.197571949935e+01, -4.742502251260e+01, 3.628013140596e+01, -1.423003075866e+01, 3.273090240144e+00, -4.557928912260e-01, 3.773588347458e-02, -1.707904867106e-03, 3.251203344615e-05};
	// Janev 1987 3.1.10
	// p + H(1s) --> H(2s) + p
	std::vector<double> aSigma_2s = {-1.327325087764e+04, 1.317576614520e+04, -5.683932157858e+03, 1.386309780149e+03, -2.089794561307e+02, 1.992976245274e+01, -1.173800576157e+00, 3.902422810767e-02, -5.606240339932e-04};

	// Contribution from ground -> 2p orbital
	double sigma_2p;
	if ( CoMEnergy < minimumEnergySigma_2p ) {
		sigma_2p = 0;
	} else {
		sigma_2p = EvaluateJanevCrossSectionFit( aSigma_2p, CoMEnergy );
	}

	// Contribution from ground -> 2s orbital
	double sigma_2s;
	if ( CoMEnergy < minimumEnergySigma_2s ) {
		sigma_2s = 0;
	} else {
		sigma_2s = EvaluateJanevCrossSectionFit( aSigma_2s, CoMEnergy );
	}

	return sigma_1s + sigma_2p + sigma_2s;
}

// Energy in electron volts, returns cross section in cm^2
double radiativeRecombinationCrossSection( double Energy )
{
	// From https://iopscience-iop-org.proxy-um.researchport.umd.edu/article/10.1088/1402-4896/ab060a
	// Igor A Kotelnikov and Alexander I Milstein 2019 Phys. Scr. 94 055403
	// Equation 9
	// H+ + e --> H + hν

	int Z = 1;
	double IonizationEnergy = 13.59844; // eV
	double J_Z = ::pow( Z, 2 ) * IonizationEnergy;
	double eta = ::sqrt( J_Z / Energy );
	double sigma = ::pow( 2, 8 ) * ::pow( M_PI * BohrRadius, 2 ) / 3 * ::pow( eta, 6 ) * ::exp( - 4 * eta * atan( 1 / eta ) ) / ( ( 1 - ::exp( -2 * M_PI * eta ) ) * ::pow( ::pow( eta, 2 ) + 1, 2 ) ) * ::pow( FineStructureConstant, 3 );

	sigma *= 1e4; // convert to cm^2
	return sigma;
}

// double electronHydrogenExcitationN2CrossSection( double Te )
// {
// 	// Minimum energy of cross section in eV
// 	const double excitationEnergy = 10.2;
// 	const double cutoffEnergy1 = 11.56;
// 	const double cutoffEnergy2 = 12.23;
//
// 	// Contribution from ground state
// 	// Janev 1993, ATOMIC AND PLASMA-MATERIAL INTERACTION DATA FOR FUSION, Volume 4
// 	// Equation 1.1.3
// 	// e + H(1s) --> e + H*(n=2)
// 	// Accuracy is 10% or better
// 	std::vector<double> fittingParamA = { 1.4182, -20.877, 49.735, -46.249, 17.442, 4.4979 };
//
// 	double sigma;
// 	if ( Te < excitationEnergy ) {
// 		sigma = 0;
// 	}
// 	else if ( Te < cutoffEnergy1 ) {
// 		sigma = 1e-16 * ( 0.255 + 0.1865 * ( Te - excitationEnergy ) );
// 	}
// 	else if ( Te < cutoffEnergy2 ) {
// 		sigma = 5.025e-17;
// 	}
// 	else {
// 		double sum = 0.0;
// 		double XEnergy = Te / excitationEnergy;
// 		for ( size_t n = 0; n < fittingParamA.size() - 1; n++ ) {
// 	      sum += fittingParamA.at( n ) / ::pow( XEnergy, n - 1 );
// 	   }
// 		sigma = 5.984e-16 / Te * ( sum + fittingParamA.back() * ::log( XEnergy ) );
// 	}
// 	return sigma;
// }
//
// double protonHydrogenExcitationN2CrossSection( double Ti )
// {
// 	// Minimum energy of cross section in keV
// 	const double minimumEnergySigma = 0.6;
// 	double TiKEV = Ti / 1000;
//
// 	// Contribution from ground state
// 	// Janev 1993, ATOMIC AND PLASMA-MATERIAL INTERACTION DATA FOR FUSION, Volume 4
// 	// Equation 2.2.1
// 	// H+ + H(1s) --> H+ + H+ + e
// 	// Accuracy is 100% or better
// 	const double A1 = 34.433;
// 	const double A2 = 44.057;
// 	const double A3 = 0.56870;
// 	const double A4 = 8.5476;
// 	const double A5 = 7.8501;
// 	const double A6 = -9.2217;
// 	const double A7 = 1.8020e-2;
// 	const double A8 = 1.6931;
// 	const double A9 = 1.9422e-3;
// 	const double A10 = 2.9068;
//
// 	double sigma;
// 	if ( TiKEV < minimumEnergySigma ) {
// 		sigma = 0;
// 	}
// 	else {
// 		// Energy is in units of keV
// 		sigma = 1e-16 * A1 * ( ::exp( -A2 / TiKEV ) * ::log( 1 + A3 * TiKEV ) / TiKEV + A4 * ::exp( -A5 * TiKEV ) / ( ::pow( TiKEV, A6 ) ) + A7 * ::exp( -A8 / TiKEV ) / ( 1 + A9 * ::pow( TiKEV, A10 ) ) );
// 	}
// 	return sigma;
// }

/*
 * Use the above functions to set steady-state neutral density/source
 *
 */
void MirrorPlasma::ComputeSteadyStateNeutrals()
{
	// Do not recalculate if we're in fixed-neutral-density mode
	if ( FixedNeutralDensity )
		return;

	// Calculate the Ionization Rate of cold neutrals from proton and electron impact:
	double IonizationRate =
		neutralsRateCoefficientCold( protonImpactIonization, *this ) * IonDensity * ReferenceDensity +
		neutralsRateCoefficientCold( electronImpactIonization, *this ) * ElectronDensity * ReferenceDensity;

#if defined( DEBUG ) && defined( ATOMIC_PHYSICS_DEBUG )
	std::cerr << "Current Ionization Rate is " << IonizationRate << std::endl;
#endif

	// Assume mix of neutrals is such that the particle densities are maintained so we just need to produce enough electrons from the source gas to balance the losses
	double ElectronLossRate = ParallelElectronParticleLoss() + ClassicalElectronParticleLosses();
	// Steady State requires
	//		Losses = n_N * n_e * IonizationRateCoefficient * Volume
	// so n_N = (Losses/Volume) / ( n_e * IonizationRateCoefficient )
	//		    = (Losses/Volume) / IonizationRate
	NeutralDensity = ElectronLossRate / ( IonizationRate );
	// Normalize neutral density
	NeutralDensity = NeutralDensity / ReferenceDensity;
	NeutralSource = ElectronLossRate;
}

// Number of ions lost as neutrals per second per unit volume
// requires computed neutral density
double MirrorPlasma::CXLossRate() const
{
	if ( !pVacuumConfig->IncludeCXLosses )
		return 0.0;

	double CXRateCoefficient = neutralsRateCoefficientCold( HydrogenChargeExchange, *this );

#if defined( DEBUG ) && defined( ATOMIC_PHYSICS_DEBUG )
	std::cerr << "Current CX Loss Rate is " << CXRateCoefficient * ( NeutralDensity * ReferenceDensity * IonDensity * ReferenceDensity ) << " particles/s"<<std::endl;
#endif

	return CXRateCoefficient * ( NeutralDensity * ReferenceDensity * IonDensity * ReferenceDensity );
}
