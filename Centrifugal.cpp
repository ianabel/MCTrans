
#include "Centrifugal.hpp"
#include "AlphaHeating.hpp"

#include <cmath>
#include <boost/math/tools/roots.hpp>

// e in Gaussian Units is e_SI / sqrt( 4 * pi * epsilon_0 )
// precomuting this to e_SI / (sqrt(4*pi))*(sqrt(epsilon_0)) gives:
constexpr double GaussianElectronCharge = ElectronCharge / ( 3.54490770181103205458 * .00000297559873134802 );

// For our purposes c_s is defined such that
// c_s^2 = Z T_e / m_i
double SoundSpeed( Plasma const& plasma )
{
	double cs = ::sqrt( plasma.ZIon*plasma.ElectronTemperature*ReferenceTemperature / ( plasma.Mu * ProtonMass ) );
	return cs;
}

double LogLambdaElectron( Plasma const& plasma )
{
	return 17.;
	// Convert to NRL Formulary Units
	double neNRL = plasma.ElectronDensity * 1e14; // We use 10^20 / m^3 = 10^14 / cm^3
	double TeNRL = plasma.ElectronTemperature * 1000; // We normalize to 1 keV they use 1 ev
	// Convert log( ne^.5 Te^-1.25 ) = .5*log(ne) - 1.25 * log(Te) 
	double logTe = ::log( TeNRL );
	return 23.5 - 0.5*::log( neNRL ) + 1.25*logTe - ::sqrt( 1e-5 + ( logTe - 2 ) * ( logTe - 2 ) / 16 );
}

// tau_ee 
double ElectronCollisionTime( Plasma const & plasma )
{
	return ::sqrt( ElectronMass * plasma.ElectronTemperature * ReferenceTemperature / 2.0 ) * ( plasma.ElectronTemperature * ReferenceTemperature ) / 
	          ( plasma.ElectronDensity * ReferenceDensity * ::pow( GaussianElectronCharge, 4 ) * LogLambdaElectron( plasma ) );
}

// See notes.
// This expression is the leading order in O(1/M^2)
double CentrifugalPotential( Plasma const& plasma, Configuration const& conf )
{
	return 0.5 * ( 1.0 - 1.0/conf.MirrorRatio ) * conf.MachNumber * conf.MachNumber / ( plasma.ZIon + plasma.TiTe );
}

// Taken from Pastukhov's paper
double ParallelElectronHeatLoss( Plasma const& plasma, Configuration const& conf )
{
	double R = conf.MirrorRatio;
	double tau_ee = ElectronCollisionTime( plasma );
	double x0 = CentrifugalPotential( plasma, conf );
	// The integral in Pastukhov is, with x = T/e Phi
	// evaluates to 1 + sqrt(pi/x) * exp(x) * erfc( sqrt(x) )
	double PastukhovIntegral = 1.0 + ::sqrt( M_PI / x0 ) * ::exp( x0 ) * std::erfc( ::sqrt( x0 ) );
	return ( M_2_SQRTPI / tau_ee ) * ( 2.0 * R / ( 2.0 * R + 1 ) ) * ( 1.0 / ::log( 4.0 * R + 2.0 ) ) * ( 2./3. + ( 1./x0 ) * PastukhovIntegral ) * ::exp( - x0 );
}


double LogLambdaIon( Plasma const& plasma )
{
	return 18.0;
	double TiNRL = plasma.ElectronTemperature * plasma.TiTe * 1000;
	double niNRL = plasma.IonDensity;
	return 23.0 - 0.5 * ::log( niNRL ) - 1.5 * ::log( plasma.ZIon * plasma.ZIon / TiNRL );
}

double IonCollisionTime( Plasma const& plasma )
{
	double PiThreeHalves = ::pow( M_PI, 1.5 ); // pi^(3/2)
	double TiThreeHalves = ::pow( plasma.ElectronTemperature * plasma.TiTe * ReferenceTemperature, 1.5 );
	return 12 * ::sqrt( plasma.Mu * ProtonMass ) * PiThreeHalves * TiThreeHalves * VacuumPermittivity * VacuumPermittivity / ( plasma.ElectronDensity * ReferenceDensity * ::pow( plasma.ZIon * ElectronCharge, 4 ) * LogLambdaIon( plasma ) );
}

// Pastukhov-style ion particle loss 
// ion potential energy drop from midplane to throat is 
// same as the potential drop, and has the SAME sign.
// So the electric potential is electron-confining, and the total PE (centrifugal + electrostatic) is ion-confining.
double ParallelIonParticleLoss( Plasma const& plasma, Configuration const& conf ) {
	return 1.;
}

double ParallelKineticEnergyLoss( Plasma const& plasma, Configuration const& conf )
{
	double R = conf.MirrorRatio;
	double tau_ii = IonCollisionTime( plasma );

	// Centrifugal potential is normalised to e/T_e
	double x0 = CentrifugalPotential( plasma, conf ) * plasma.ZIon / plasma.TiTe;

	double PastukhovIntegral = 1 + ::sqrt( M_PI / x0 ) * ::exp( x0 ) * std::erfc( ::sqrt( x0 ) );

	// This is a bad estimate, as some deconfining will occur to ensure quasineutrality.
	double ParticleLossRate = ( M_2_SQRTPI / tau_ii ) * ( 2.0 * R / ( 2.0 * R + 1 ) ) * ( 1.0 / ::log( 4.0 * R + 2.0 ) ) * ( 2./3. + ( 1./x0 ) * PastukhovIntegral ) * ::exp( - x0 );

	double EndExpansionRatio = 1.0; // If particles are hitting the wall at a radius larger than they are confined at, this increases momentum loss per particle 
	double KineticEnergyPerParticle = 0.5 * plasma.ElectronTemperature * ReferenceTemperature * ( EndExpansionRatio * EndExpansionRatio ) *  conf.MachNumber * conf.MachNumber;
	// Convert to W/m^3 for consistency across loss rates
	return ParticleLossRate * plasma.ElectronDensity * ReferenceDensity * KineticEnergyPerParticle;
};

double IonCyclotronFrequency( Plasma const& plasma, Configuration const& conf )
{
	double MagneticField = conf.CentralCellFieldStrength;
	return plasma.ZIon * ElectronCharge * MagneticField / ( plasma.Mu * ProtonMass );
}


double ClassicalHeatLossRate( Plasma const& plasma, Configuration const& conf )
{
	double omega_ci = IonCyclotronFrequency( plasma, conf );
	double kappa_perp = 2.0 * plasma.ElectronTemperature * ReferenceTemperature * plasma.TiTe / ( plasma.Mu * ProtonMass * omega_ci * omega_ci * IonCollisionTime( plasma ) );
	double L_T = ( conf.PlasmaColumnWidth / 2.0 );
	return kappa_perp /( L_T * L_T );
}

// Braginskii eta_1
double ClassicalViscosity( Plasma const& plasma, Configuration const& conf )
{
	double omega_ci = IonCyclotronFrequency( plasma, conf );
	return ( 3.0 / 10.0 ) * ( plasma.ZIon * plasma.ElectronDensity * ReferenceDensity * plasma.TiTe * plasma.ElectronTemperature * ReferenceTemperature ) / ( omega_ci * omega_ci * IonCollisionTime( plasma ) );
}

// Viscous heating = eta * u^2 / L_u^2
double ViscousHeating( Plasma const& plasma, Configuration const& conf )
{
	double L_u = ( conf.PlasmaColumnWidth / 2.0 );
	double Velocity = conf.MachNumber * SoundSpeed( plasma );
	double VelocityShear = Velocity / L_u;

	return ClassicalViscosity( plasma, conf ) * VelocityShear * VelocityShear;
}

double AlfvenMachNumber( Plasma const& plasma, Configuration const& conf )
{
	// v_A^2 = B^2 / mu_0 * n_i * m_i
	// c_s^2 = T_e / m_i
	double IonDensity = plasma.ElectronDensity * ReferenceDensity / plasma.ZIon;
	double MassDensity = IonDensity * ( plasma.Mu * ProtonMass );
	double AlfvenSpeed = conf.CentralCellFieldStrength / ::sqrt( PermeabilityOfFreeSpace * MassDensity );

	return conf.MachNumber * ( SoundSpeed( plasma ) / AlfvenSpeed );
}

double PlasmaVolume( Configuration const& conf )
{
	return M_PI * ( conf.PlasmaColumnWidth + 2 * conf.AxialGapDistance ) * conf.PlasmaColumnWidth * conf.PlasmaLength;
}

double WallArea( Configuration const& conf )
{
	return 2.0 * M_PI * ( conf.WallRadius * conf.WallRadius ) * conf.PlasmaLength; 
}

double BremsstrahlungLosses( Plasma const& plasma )
{
	// NRL formulary with reference values factored out
	return 169 * ::sqrt( 1000 * plasma.ElectronTemperature ) * plasma.Zeff * plasma.Zeff * plasma.ElectronDensity * plasma.ElectronDensity;
}

// Approximate tau_E^-1 by the sum of the heat loss rates

double NeutralLossRate( Plasma const& plasma, Configuration const& conf )
{
	return 0.0;
	double NeutralDensity;
	// Plasma velocity
	double u = conf.MachNumber * SoundSpeed( plasma );

	if ( plasma.NeutralSource > 0 ) {
		// Work out steady-state neutral density
		// N.B. CX causes momentum loss but no increase in neutral density
		// Ionization rate per neutral
		double IonizationRate = IonizationCrossSection * u;
	} else {
		NeutralDensity = plasma.NeutralDensity;
	}

	// Mean free path of plasma particle in cloud of neutrals
	double cx_mfp = NeutralDensity * CXCrossSection; 
	// Loss rate per particle due to CX
	double cx_rate = u / cx_mfp;
	double IonDensity = plasma.ElectronDensity / plasma.ZIon;
	double TotalParticleLossRate = IonDensity * PlasmaVolume( conf )*cx_rate;

	return TotalParticleLossRate;
}

double TauEInverse( Plasma const& plasma, Configuration const& conf )
{
	return ( ClassicalHeatLossRate( plasma, conf ) + ParallelElectronHeatLoss( plasma, conf ) + NeutralLossRate( plasma, conf ) );
}

// add Bremsstrahlung radiative losses
double TotalHeatLoss( Plasma const& plasma, Configuration const& conf )
{

	double IonDensity = plasma.ElectronDensity / plasma.ZIon;
	double IonTemperature = plasma.ElectronTemperature * plasma.TiTe;

	double StoredEnergy = 1.5 * ( plasma.ElectronDensity * plasma.ElectronTemperature + IonDensity * IonTemperature ) * ReferenceDensity * ReferenceTemperature;
	
	return StoredEnergy * TauEInverse( plasma, conf ) + BremsstrahlungLosses( plasma );
}


// Steady State Mach Number
void SetMachNumber( Plasma const &plasma, Configuration &conf )
{
	// NB This uses power densities in W/m^3
	auto PowerBalance = [ &plasma, &conf ]( double M ) {
		Configuration test_conf = conf;
		test_conf.MachNumber = M;
		double HeatLoss = TotalHeatLoss( plasma, test_conf );
		double Heating = ViscousHeating( plasma, test_conf ) + test_conf.AuxiliaryHeating*1e6 / PlasmaVolume( test_conf );
		if ( conf.IncludeAlphaHeating ) {
			/*
			 * If the slowing-down time is longer than the energy confinement time
			 * reduce the alpha heating by that amount
			 */
			double AlphaHeatingFraction = 1. / ( SlowingDownTime( plasma ) * TauEInverse( plasma, conf ) );
			// if ( AlphaHeatingFraction > 1.0 )
			AlphaHeatingFraction = 1.0;
			// 1e6 as we are using W/m^3 and the formulary was in MW/m^3
			Heating += AlphaHeating( plasma, test_conf ) * 1e6 * AlphaHeatingFraction;

		}
		return Heating - HeatLoss;
	};

	boost::uintmax_t iters = 1000;
	boost::math::tools::eps_tolerance<double> tol( 11 ); // only bother getting part in 1024 accuracy
	double InitialMach = 1.0; // Usually M > 4 for these solutions
	double Factor = 2.0;
	bool rising = true; // Confinement gets uniformly better for increasing M, and Viscous heating increases with M
	auto [ M_lower, M_upper ] = boost::math::tools::bracket_and_solve_root( PowerBalance, InitialMach, Factor, rising, tol, iters );
	conf.MachNumber = ( M_lower + M_upper )/2.0;
}

void SetMachFromVoltage( Plasma const& plasma, Configuration &conf )
{
	// u = E x B / B^2 
	// M = u/c_s ~ (V/aB)/cs
	double cs = SoundSpeed( plasma );
	double Mach = conf.Voltage / ( conf.PlasmaColumnWidth * conf.CentralCellFieldStrength * cs );
	conf.MachNumber = Mach;
}

void SetTemperature( Plasma &plasma, Configuration &conf )
{
	// NB This uses power densities in W/m^3
	auto PowerBalance = [ &plasma, &conf ]( double Te ) {
		Plasma test_plasma = plasma;
		test_plasma.ElectronTemperature = Te;
		// Update Mach Number from new T_e
		SetMachFromVoltage( test_plasma, conf );
		double HeatLoss = TotalHeatLoss( test_plasma, conf );
		double Heating = ViscousHeating( test_plasma, conf ) + conf.AuxiliaryHeating*1e6 / PlasmaVolume( conf );
		if ( conf.IncludeAlphaHeating ) {
			/*
			 * If the slowing-down time is longer than the energy confinement time
			 * reduce the alpha heating by that amount
			 */
			double AlphaHeatingFraction = 1. / ( SlowingDownTime( test_plasma ) * TauEInverse( test_plasma, conf ) );
			// if ( AlphaHeatingFraction > 1.0 )
			AlphaHeatingFraction = 1.0;
			// 1e6 as we are using W/m^3 and the formulary was in MW/m^3
			Heating += AlphaHeating( test_plasma, conf ) * 1e6 * AlphaHeatingFraction;

		}
		return Heating - HeatLoss;
	};

	boost::uintmax_t iters = 1000;
	boost::math::tools::eps_tolerance<double> tol( 11 ); // only bother getting part in 1024 accuracy
	double InitialTe = 0.1; // Start cold -- 100eV
	double Factor = 2.0;
	bool rising = false; // At fixed Mach Number, heating the plasma up will increase losses
	auto [ T_lower, T_upper ] = boost::math::tools::bracket_and_solve_root( PowerBalance, InitialTe, Factor, rising, tol, iters );
	plasma.ElectronTemperature = ( T_lower + T_upper )/2.0;
	SetMachFromVoltage( plasma, conf );
}

double ElectricPotential( Plasma const &plasma, Configuration const& conf )
{
	double RotationVelocity = SoundSpeed( plasma )*conf.MachNumber;
	return conf.CentralCellFieldStrength * RotationVelocity * conf.PlasmaColumnWidth; 
}

double IonLarmorRadius( Plasma const& plasma, Configuration const& conf )
{
	double Ti = plasma.TiTe * plasma.ElectronTemperature * 1000;
	return 1.02 * ::sqrt( plasma.Mu * Ti ) / ( plasma.ZIon * conf.CentralCellFieldStrength * 10000 );
}

double Beta( Plasma const& plasma, Configuration const& conf )
{
	// From Plasma Formulary, so convert to /cm^3 , eV, and Gauss
	double ne_Formulary = plasma.ElectronDensity * 1e-6; 
	double Te_Formulary = plasma.ElectronTemperature * 1e3;
	double IonDensity = plasma.ElectronDensity * ReferenceDensity / plasma.ZIon;
	double ni_Formulary = IonDensity * 1e-6;
	double Ti_Formulary = plasma.ElectronTemperature * plasma.TiTe * 1e3;
	double MagField_Formulary = conf.CentralCellFieldStrength * 1e4;
	return 4.03e-11 * ( ne_Formulary * Te_Formulary + ni_Formulary * Ti_Formulary ) / ( MagField_Formulary * MagField_Formulary );
}

double DebyeLength( Plasma const& plasma )
{
	// l_debye^2 = epsilon_0 * (k_B T_e) / n_e * e^2
	double LambdaDebyeSquared = VacuumPermittivity * ( plasma.ElectronTemperature * ReferenceTemperature ) / ( plasma.ElectronDensity * ReferenceDensity * ElectronCharge * ElectronCharge );
	return ::sqrt( LambdaDebyeSquared );
}


double NuStar( Plasma const& plasma, Configuration const& conf )
{
	return conf.PlasmaLength / ( IonCollisionTime( plasma ) * SoundSpeed( plasma ) );
}


double KineticEnergy( Plasma const& plasma, Configuration const& conf )
{
	double IonMass = ProtonMass * plasma.Mu;
	double u = conf.MachNumber * SoundSpeed( plasma );
	double IonDensity = plasma.ElectronDensity / plasma.ZIon;
	return .5 * IonMass * IonDensity * ReferenceDensity * u * u * PlasmaVolume( conf );
}

double ThermalEnergy( Plasma const& plasma, Configuration const& conf )
{
	double IonDensity = plasma.ElectronDensity / plasma.ZIon;
	double IonTemperature = plasma.ElectronTemperature * plasma.TiTe;
	return 1.5 * ( plasma.ElectronDensity * plasma.ElectronTemperature + IonDensity * IonTemperature ) * ReferenceDensity * ReferenceTemperature * PlasmaVolume( conf ); 
}


