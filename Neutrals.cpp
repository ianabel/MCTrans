
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


Species Electron{ .type = Species::Electron, .Charge = -1, .Mass = ElectronMass,  .Name = "Electron" };
Species Proton{ .type = Species::Ion, .Charge = 1, .Mass = ProtonMass, .Name = "Proton"};
Species Deuteron{ .type = Species::Ion, .Charge = 1, .Mass = 1.999*ProtonMass, .Name = "Deuteron" };
Species NeutralHydrogen{ .type = Species::Neutral, .Charge = 0, .Mass = ProtonMass + ElectronMass, .Name = "Neutral Hydrogen"};

// Cross Section objects for integration
CrossSection protonImpactIonization( protonImpactIonizationCrossSection, 200, 1e6, Proton, NeutralHydrogen );
CrossSection HydrogenChargeExchange( HydrogenChargeExchangeCrossSection, 0.1, 1e6, Proton, NeutralHydrogen );
CrossSection electronImpactIonization( electronImpactIonizationCrossSection, 13.6, 1e6, Electron, NeutralHydrogen );

/* We don't actually use this....
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

	constexpr double tolerance = 1e-6;
	constexpr unsigned MaxDepth = 10;
	double HotRateCoeff = 4.0 / ( ::sqrt(2 * M_PI * sigma.Particle.Mass * temperature) * temperature )
	         * boost::math::quadrature::gauss_kronrod<double, 15>::integrate( integrand, sigma.MinEnergy, sigma.MaxEnergy, MaxDepth, tolerance );

#if defined( DEBUG ) && defined( ATOMIC_PHYSICS_DEBUG )
	std::cerr << "Computing a hot rate coefficient at T = " << plasma.ElectronTemperature/1000 << " eV and M = " << plasma.MachNumber << " gave <sigma v> = " << HotRateCoeff  << std::endl;
#endif

	return HotRateCoeff;

}
*/

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

	int delta_ns = 0;
	if ( sigma.Particle.Name == sigma.Target.Name ){
		delta_ns = 1;
	}

	double thermalSpeed = ::sqrt( 2.0 * temperature / sigma.Particle.Mass );
	double thermalMachNumber = plasma.MachNumber * plasma.SoundSpeed() / thermalSpeed;

	// Needs Energy in electronVolts
	auto integrand = [&]( double Energy ) {
		double velocity = ::sqrt( 2.0 * Energy * ElectronCharge / sigma.Particle.Mass );
		double Jacobian = ElectronCharge / ( sigma.Particle.Mass * velocity ); // The integral is over Energy, which is in units of electronvolts, so transform the integrand back to eV, including change of variables from dvelocity to dE
		double sigmaM2 = sigma( Energy ) * 1e-4; // sigma is in cm^2, we need m^2
		double MV1 = thermalMachNumber - velocity / thermalSpeed;
		double MV2 = thermalMachNumber + velocity / thermalSpeed;
		return velocity * velocity * sigmaM2 * ( ::exp( - MV1*MV1 ) - ::exp( -MV2*MV2 ) ) * Jacobian;
	};

	constexpr double tolerance = 1e-6;
	constexpr unsigned MaxDepth = 10;


	double min_sqrt = 4;
	// Velocity such that Exp( -(M-v)^2 ) is negligibly small;
	double min_velocity;
	if ( thermalMachNumber <= min_sqrt )
		min_velocity = 0;
	else 
		min_velocity = thermalMachNumber - min_sqrt;
	double max_velocity;
	max_velocity = thermalMachNumber + min_sqrt;

	double minEnergy = std::max( sigma.MinEnergy, min_velocity*min_velocity * temperature / ElectronCharge );
	double maxEnergy = std::min( sigma.MaxEnergy, max_velocity * max_velocity * temperature / ElectronCharge );

	double ColdRateCoeff = 1 / ( thermalMachNumber * thermalSpeed * thermalSpeed * ::sqrt(M_PI) * ( 1 + delta_ns ) )
	        * boost::math::quadrature::gauss_kronrod<double, 15>::integrate( integrand, minEnergy, maxEnergy, MaxDepth, tolerance );

#if defined( DEBUG ) && defined( ATOMIC_PHYSICS_DEBUG )
	std::cerr << "Computing a cold rate coefficient at T = " << 1000*temperature/ReferenceTemperature  << " eV and M = " << plasma.MachNumber << " gave <sigma v> = " << ColdRateCoeff  << std::endl;
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
	return density_1 * density_2 * rateCoefficient * pMirrorPlasma->PlasmaVolume();
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

double electronImpactIonizationCrossSection( double Energy )
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
	if ( Energy < minimumEnergySigma ) {
		return 0.0;
	}
	else {
		double sum = 0.0;
		double x = 1.0 - ionizationEnergy / Energy;
		if ( x <= 0 )
			return 0.0;
		for ( size_t n = 0; n < fittingParamB.size(); n++ ) {
	      sum += fittingParamB.at( n ) * ::pow( x, n+1 );
	   }
		sigma = ( 1.0e-13 / ( ionizationEnergy * Energy ) ) * ( fittingParamA * ::log( Energy / ionizationEnergy ) + sum );
		return sigma;
	}
}

// Energy in electron volts, returns cross section in cm^2
double protonImpactIonizationCrossSection( double Energy )
{
	// Minimum energy of cross section in keV
	const double minimumEnergySigma = 0.5;
	// Convert to keV
	double EnergyKEV = Energy / 1000;

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
	if ( EnergyKEV < minimumEnergySigma ) {
		sigma = 0;
	}
	else {
		// Energy is in units of keV
		sigma = 1e-16 * A1 * ( ::exp( -A2 / EnergyKEV ) * ::log( 1 + A3 * EnergyKEV ) / EnergyKEV + A4 * ::exp( -A5 * EnergyKEV ) / ( ::pow( EnergyKEV, A6 ) + A7 * ::pow( EnergyKEV, A8 ) ) );
	}
	return sigma;
}

double HydrogenChargeExchangeCrossSection( double Energy )
{
	// Minimum energy of cross section in eV
	constexpr double minimumEnergySigma_n1 = 0.12;
	// constexpr double minimumEnergySigma_n2 = 10;
	// constexpr double minimumEnergySigma_n3 = 10;

	// Contribution from ground -> ground state
	// Janev 1993 2.3.1
	// p + H(n=1) --> H + p
	double sigma_n1;
	if ( Energy < minimumEnergySigma_n1 ) {
		sigma_n1 = 0;
	} else {
		double EnergyKEV = Energy / 1000;
		sigma_n1 = 1e-16 * 3.2345 * ::log( 235.88 / EnergyKEV + 2.3713 ) / ( 1 + 0.038371 * EnergyKEV + 3.8068e-6 * ::pow( EnergyKEV, 3.5 ) + 1.1832e-10 * ::pow( EnergyKEV, 5.4 ) );
	}

	/*
	// Contribution from n=2 orbital
	double sigma_n2;
	if ( Energy < minimumEnergySigma_n2 ) {
		sigma_n2 = 0;
	} else {
		double EnergyKEV = Energy / 1000;
		int n = 2;
		double EnergyTilde = EnergyKEV * ::pow( n, 2 );
		sigma_n2 = 1e-16 * 0.92750 * ::log( 6.5040e3 / EnergyTilde + 20.699 ) / ( 1 + 1.3405e-2 * EnergyTilde + 3.0842e-6 * ::pow( EnergyTilde, 3.5 ) + 1.1832e-10 * ::pow( EnergyTilde, 5.4 ) );
	}

	// Contribution from n=3 orbital
	double sigma_n3;
	if ( Energy < minimumEnergySigma_n3 ) {
		sigma_n3 = 0;
	} else {
		double EnergyKEV = Energy / 1000;
		int n = 3;
		double EnergyTilde = EnergyKEV * ::pow( n, 2 );
		sigma_n3 = 1e-16 * 0.37271 * ::log( 2.7645e6 / EnergyTilde + 1.4857e3 ) / ( 1 + 1.5720e-3 * EnergyTilde + 3.0842e-6 * ::pow( EnergyTilde, 3.5 ) + 1.1832e-10 * ::pow( EnergyTilde, 5.4 ) );
	}
	*/

	// IGA: I believe with cold neutrals the contribution from higher orbitals is negligible
	// but leaving as-is for the moment as a pessimistic assumption
	// NRS: I changed the cross section to only include the ground state
	return sigma_n1;
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
	if ( !IncludeCXLosses )
		return 0.0;

	double CXRateCoefficient = neutralsRateCoefficientCold( HydrogenChargeExchange, *this );

	// Really rate per unit plasma volume
	double CXRate = CXRateCoefficient * ( NeutralDensity * ReferenceDensity * IonDensity * ReferenceDensity );

#if defined( DEBUG ) && defined( ATOMIC_PHYSICS_DEBUG )
	std::cerr << "Current CX Loss Rate is " << CXRate << " particles/s/m^3"<<std::endl;
	std::cerr << " CX Frequency is " << CXRateCoefficient * NeutralDensity * ReferenceDensity << " /s " << std::endl;
#endif

	return CXRate;
}
