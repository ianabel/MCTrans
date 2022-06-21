

#include "MirrorPlasma.hpp"
#include "PlasmaPhysics.hpp"
#include <cmath>

/*
 * We assume pure mirror confinement for Alpha particles.
 * Hence if the mirror ration is R, all particles with
 * v_perp / v < sqrt( 1/R ) will be lost along the field line
 */

// Born at 3.52 MeV
constexpr double AlphaBirthEnergy = 3.52e6 * ElectronCharge;

double MirrorPlasma::PromptAlphaLossFraction() const
{
	// Transforming to pitch angle and integrating, 
	// the formula is
	// Loss Fraction = 1 - Sqrt(1 - 1/R)
	return 1.0 - ::sqrt( 1.0 - 1.0 / MirrorRatio );
}

// Just use classical formula 
// ln Lambda = ln ( Debye Length / Closest Approach )
double MirrorPlasma::LogLambdaAlphaElectron() const
{
	double AlphaElectronReducedMass = 1. / ( 1./ElectronMass + 1./AlphaMass );
	double AlphaEnergy = 3.52 * 1e6 * ElectronCharge;
	double AlphaVelocity = ::sqrt( 2.0 * AlphaEnergy / AlphaMass );
	// Assume fast alphas on motionless electrons
	double RelativeVelocity = AlphaVelocity;
	double ClosestApproach = Z_Alpha * ElectronCharge * ElectronCharge / ( AlphaElectronReducedMass * RelativeVelocity );
	return ::log( DebyeLength() / ClosestApproach );
}

// Tau_SD for alpha particles
double MirrorPlasma::SlowingDownTime() const
{
	// From Helander & Sigmar, boxed equation after (3.50)
	// 3 * (2 pi)^(3/2) * Epsilon_0 ^2 * mAlpha * ( ElectronTemperature )^(3/2) /
	//   Z_alpha^2 e^4 m_e^(1/2) * n_e * log Lambda
	return 3.0 * ::pow( 2.0 * M_PI, 1.5 ) * ( VacuumPermittivity  * VacuumPermittivity * AlphaMass * ::pow( ElectronTemperature * ReferenceTemperature, 1.5 ) ) /
			( Z_Alpha * Z_Alpha * ElectronCharge * ElectronCharge * ElectronCharge * ElectronCharge * ::sqrt( ElectronMass ) * ElectronDensity * ReferenceDensity * LogLambdaAlphaElectron() );
}

// We assume alphas are roughly lost as if they are a species
// at T_alpha = m_alpha v_crit^2 / 2
// experiencing pitch-angle scattering at a rate 1/tau_s
double MirrorPlasma::AlphaParallelLossRate() const
{
	double tau_SD = SlowingDownTime();

	double R = MirrorRatio;

	// In Units of T_e
	double CriticalEnergy = ::pow( IonSpecies.Charge * 3.0 * ::sqrt( M_PI ) * ElectronMass / ( IonSpecies.Mass * ProtonMass ), 2./3. ) * ( AlphaMass / ElectronMass ); 


	double x0 = CentrifugalPotential() * Z_Alpha / CriticalEnergy;
	double PastukhovIntegral = 1 + ::sqrt( M_PI / x0 ) * ::exp( x0 ) * std::erfc( ::sqrt( x0 ) );

	// This is a bad estimate (relies on R>>1)
	//
	double ParticleLossRate = ( M_2_SQRTPI / tau_SD ) * ( 2.0 * R / ( 2.0 * R + 1 ) ) * ( 1.0 / ::log( 4.0 * R + 2.0 ) ) * ( 2./3. + ( 1./x0 ) * PastukhovIntegral ) * ::exp( - x0 );
	return ParticleLossRate;
}

double MirrorPlasma::AlphaHeating() const
{
	// double tau_SD = SlowingDownTime( plasma );
	double LossFraction = PromptAlphaLossFraction();
	double GrossHeatingRate = FusionAlphaPowerDensity();

	// Currently just assume all particles not lost promptly are confined.
	return GrossHeatingRate * ( 1.0 - LossFraction );
}

double MirrorPlasma::AlphaPromptLosses() const
{
	double LossFraction = PromptAlphaLossFraction();
	double GrossHeatingRate = FusionAlphaPowerDensity();

	// Currently just assume all particles not lost promptly are confined.
	return GrossHeatingRate * LossFraction;
}

// I.e. Parallel momentum loss rate.
double MirrorPlasma::PromptAlphaThrust() const
{
	// Integrate m_alpha v_|| * isotropic alpha birth distribution
	// Lost momentum = Alpha Birth Rate * Alpha Birth Momentum / Mirror Ratio
	double AlphaMomentum = ::sqrt( 2.0 * AlphaMass * AlphaBirthEnergy );
	double AlphasPerSecond = AlphaProductionRate() * PlasmaVolume();
	return AlphasPerSecond * AlphaMomentum / MirrorRatio;
}
