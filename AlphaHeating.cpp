

#include "AlphaHeating.hpp"
#include <cmath>


// CODATA value 6.644657230(82)×10−27 kg
double AlphaMass = 6.644657230e-27;
double Z_Alpha = 2.0;

/*
 * We assume pure mirror confinement for Alpha particles.
 * Hence if the mirror ration is R, all particles with
 * v_perp / v < sqrt( 1/R ) will be lost along the field line
 */


double PromptAlphaLossFraction( Plasma const& plasma, Configuration const& conf )
{
	// Transforming to pitch angle and integrating, 
	// the formula is
	// Loss Fraction = 1 - Sqrt(1 - 1/R)
	return 1.0 - ::sqrt( 1.0 - 1.0 / conf.MirrorRatio );
}

// Just use classical formula 
// ln Lambda = ln ( Debye Length / Closest Approach )
double LogLambdaAlphaElectron( Plasma const& plasma )
{
	double AlphaElectronReducedMass = 1. / ( 1./ElectronMass + 1./AlphaMass );
	double AlphaEnergy = 3.52 * 1e6 * ElectronCharge;
	double AlphaVelocity = ::sqrt( 2.0 * AlphaEnergy / AlphaMass );
	// Assume fast alphas on motionless electrons
	double RelativeVelocity = AlphaVelocity;
	double ClosestApproach = Z_Alpha * ElectronCharge * ElectronCharge / ( AlphaElectronReducedMass * RelativeVelocity );
	return ::log( DebyeLength( plasma ) / ClosestApproach );
}

double SlowingDownTime( Plasma const& plasma )
{
	// From Helander & Sigmar, boxed equation after (3.50)
	// 3 * (2 pi)^(3/2) * Epsilon_0 ^2 * mAlpha * ( ElectronTemperature )^(3/2) /
	//   Z_alpha^2 e^4 m_e^(1/2) * n_e * log Lambda
	return 3.0 * ::pow( 2.0 * M_PI, 1.5 ) * ( VacuumPermittivity  * VacuumPermittivity * AlphaMass * ::pow( plasma.ElectronTemperature * ReferenceTemperature, 1.5 ) ) /
			( Z_Alpha * Z_Alpha * ElectronCharge * ElectronCharge * ElectronCharge * ElectronCharge * ::sqrt( ElectronMass ) * plasma.ElectronDensity * ReferenceDensity * LogLambdaAlphaElectron( plasma ) );
}

// We assume alphas are roughly lost as if they are a species
// at T_alpha = m_alpha v_crit^2 / 2
// experiencing pitch-angle scattering at a rate 1/tau_s
double AlphaParallelLossRate( Plasma const& plasma, Configuration const& conf )
{
	double tau_SD = SlowingDownTime( plasma );

	double R = conf.MirrorRatio;

	// In Units of T_e
	double CriticalEnergy = ::pow( plasma.ZIon * 3.0 * ::sqrt( M_PI ) * ElectronMass / ( plasma.Mu * ProtonMass ), 2./3. ) * ( AlphaMass / ElectronMass ); 


	double x0 = CentrifugalPotential( plasma, conf ) * Z_Alpha / CriticalEnergy;
	double PastukhovIntegral = 1 + ::sqrt( M_PI / x0 ) * ::exp( x0 ) * std::erfc( ::sqrt( x0 ) );

	// This is a bad estimate (relies on R>>1)
	//
	double ParticleLossRate = ( M_2_SQRTPI / tau_SD ) * ( 2.0 * R / ( 2.0 * R + 1 ) ) * ( 1.0 / ::log( 4.0 * R + 2.0 ) ) * ( 2./3. + ( 1./x0 ) * PastukhovIntegral ) * ::exp( - x0 );
	return ParticleLossRate;
}

double AlphaHeating( Plasma const& plasma, Configuration const& conf )
{
	// double tau_SD = SlowingDownTime( plasma );
	double LossFraction = PromptAlphaLossFraction( plasma, conf );
	double GrossHeatingRate = FusionAlphaPowerDensity( plasma, conf );

	// Currently just assume all particles not lost promptly are confined.
	return GrossHeatingRate * ( 1.0 - LossFraction );
}

double AlphaPromptLosses( Plasma const& plasma, Configuration const& conf )
{
	double LossFraction = PromptAlphaLossFraction( plasma, conf );
	double GrossHeatingRate = FusionAlphaPowerDensity( plasma, conf );

	// Currently just assume all particles not lost promptly are confined.
	return GrossHeatingRate * LossFraction;
}


