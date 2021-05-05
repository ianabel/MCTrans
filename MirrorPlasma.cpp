
#include "MirrorPlasma.hpp"
#include <cmath>
#include <iostream>

MirrorPlasma::VacuumMirrorConfiguration::VacuumMirrorConfiguration( toml::value const& plasmaConfig )
{
	if ( plasmaConfig.count( "IonSpecies" ) != 1 )
		throw std::invalid_argument( "Fuel must be specified once in the [plasma] block" );

	std::string FuelName = plasmaConfig.at( "IonSpecies" ).as_string();

	if ( FuelName == "Hydrogen" ) {
		IonSpecies.Mass   = 1.0;
		IonSpecies.Charge = 1.0;
		IonSpecies.Name   = "Hydrogen";
		AlphaHeating = false;
	} else if ( FuelName == "Deuterium" ) {
		IonSpecies.Mass   = 2.0;
		IonSpecies.Charge = 1.0;
		IonSpecies.Name   = "Deuterium";
		AlphaHeating = false;
	} else if ( FuelName == "DT Fuel" ) {
		IonSpecies.Mass   = 2.5;
		IonSpecies.Charge = 1.0;
		IonSpecies.Name   = "Deuterium/Tritium Fuel";
		AlphaHeating = true;
	} else {
		std::string ErrorMessage = "Fuel is not a recognized plasma species";
		throw std::invalid_argument( ErrorMessage );
	}

	const auto mirrorConfig = toml::find<toml::table>( plasmaConfig, "configuration" );


	CentralCellFieldStrength =  mirrorConfig.at( "CentralCellField" ).as_floating();

	if ( mirrorConfig.count( "MirrorRatio" ) == 1 ) {
		// set Mirror Ratio directly
		MirrorRatio = mirrorConfig.at( "MirrorRatio" ).as_floating();

		if ( mirrorConfig.count( "ThroatField" ) == 1 ) {
			throw std::invalid_argument( toml::format_error( "[error] Cannot specify botth Miror Ratio and Throat Field",mirrorConfig.at( "MirrorRatio" )," mirror ratio defined here", mirrorConfig.at( "ThroatField" ), " Throat field here" ) );
		}
	} else if ( mirrorConfig.count( "ThroatField" ) == 1 ) {
		double MagFieldThroat = mirrorConfig.at( "ThroatField" ).as_floating();
		MirrorRatio = MagFieldThroat / CentralCellFieldStrength;
		if ( mirrorConfig.count( "MirrorRatio" ) == 1 ) {
			throw std::invalid_argument( toml::format_error( "[error] Cannot specify botth Miror Ratio and Throat Field",mirrorConfig.at( "MirrorRatio" )," mirror ratio defined here", mirrorConfig.at( "ThroatField" ), " Throat field here" ) );
		}
	} else {
		throw std::invalid_argument( "[error] Must specify either MirrorRatio or ThroatField" );
	}

	if ( mirrorConfig.count( "PlasmaRadiusMin" ) == 1  || mirrorConfig.count( "PlasmaRadiusMax" ) == 1 ) {
		if ( mirrorConfig.count( "PlasmaRadiusMin" ) == 0 )
		{
			throw std::invalid_argument( toml::format_error( "[error] When PlasmaRadiusMax is specified you must also set PlasmaRadiusMin",mirrorConfig.at( "PlasmaRadiusMax" )," max radius set here" ) );
		}

		if ( mirrorConfig.count( "PlasmaRadiusMax" ) == 0 )
		{
			throw std::invalid_argument( toml::format_error( "[error] When PlasmaRadiusMin is specified you must also set PlasmaRadiusMax",mirrorConfig.at( "PlasmaRadiusMin" )," min radius set here" ) );
		}
		// set Plasma Radii directly
		double PlasmaMinRadius = mirrorConfig.at( "PlasmaMinRadius" ).as_floating();
		double PlasmaMaxRadius = mirrorConfig.at( "PlasmaMaxRadius" ).as_floating();

		AxialGapDistance = PlasmaMinRadius;
		PlasmaColumnWidth = PlasmaMaxRadius - PlasmaMinRadius;


		if ( mirrorConfig.count( "AxialGapDistance" ) == 1 ) {
			throw std::invalid_argument( toml::format_error( "[error] When PlasmaRadiusMin / Max are specified you cannot set AxialGapDistance",mirrorConfig.at( "PlasmaRadiusMin" )," minimum radius set here", mirrorConfig.at( "AxialGapDistance" ), " AxialGapDistance found here" ) );
		} else if ( mirrorConfig.count( "PlasmaColumnWidth") == 1 ) {
			throw std::invalid_argument( toml::format_error( "[error] When PlasmaRadiusMin / Max are specified you cannot set AxialGapDistance",mirrorConfig.at( "PlasmaRadiusMin" )," minimum radius set here", mirrorConfig.at( "AxialGapDistance" ), " AxialGapDistance found here" ) );
		}
	} else if ( mirrorConfig.count( "AxialGapDistance" ) == 1  || mirrorConfig.count( "PlasmaColumnWidth" ) == 1 ) {
		if ( mirrorConfig.count( "AxialGapDistance" ) == 0 )
		{
			throw std::invalid_argument( toml::format_error( "[error] When PlasmaColumnWidth is specified you must also set AxialGapDistance",mirrorConfig.at( "PlasmaColumnWidth" )," width set here" ) );
		}

		if ( mirrorConfig.count( "PlasmaColumnWidth" ) == 0 )
		{
			throw std::invalid_argument( toml::format_error( "[error] When AxialGapDistance is specified you must also set PlasmaColumnWidth",mirrorConfig.at( "AxialGapDistance" )," axial gap set here" ) );
		}
		// set parameters directly

		AxialGapDistance = mirrorConfig.at( "AxialGapDistance" ).as_floating();
		PlasmaColumnWidth = mirrorConfig.at( "PlasmaColumnWidth" ).as_floating();

		if ( mirrorConfig.count( "PlasmaMinRadius" ) == 1 ) {
			throw std::invalid_argument( toml::format_error( "[error] When AxialGapDistance and PlasmaColumnWidth are specified you cannot set PlasmaMinRadius",mirrorConfig.at( "PlasmaRadiusMin" )," minimum radius set here", mirrorConfig.at( "AxialGapDistance" ), " AxialGapDistance found here" ) );
		} else if ( mirrorConfig.count( "PlasmaMaxRadius") == 1 ) {
			throw std::invalid_argument( toml::format_error( "[error] When AxialGapDistance and PlasmaColumnWidth are specified you cannot set PlasmaMaxRadius",mirrorConfig.at( "PlasmaRadiusMax" )," maximum radius set here", mirrorConfig.at( "AxialGapDistance" ), " AxialGapDistance found here" ) );
		}
	} else {
		throw std::invalid_argument( "[error] Must specify plasma annulus with either PlasmaRadiusMin / PlasmaRadiusMax or AxialGapDistance and PlasmaColumnWidth" );
	}

	if ( mirrorConfig.count( "Voltage" ) == 1 )
		ImposedVoltage = mirrorConfig.at( "Voltage" ).as_floating();
	else if ( mirrorConfig.count( "Voltage" ) == 0 )
		ImposedVoltage = 0.0;
	else 
		throw std::invalid_argument( "Imposed voltage specified more than once!" );

	WallRadius = mirrorConfig.at( "WallRadius" ).as_floating();

	PlasmaLength = mirrorConfig.at( "PlasmaLength" ).as_floating();

	if ( mirrorConfig.count( "AuxiliaryHeating" ) == 1 )
		AuxiliaryHeating = mirrorConfig.at( "AuxiliaryHeating" ).as_floating();
	else
		AuxiliaryHeating = 0.0;

	// This overrides the default for the chosen fuel
	if ( mirrorConfig.count( "IncludeAlphaHeating" ) == 1 )
		AlphaHeating = mirrorConfig.at( "IncludeAlphaHeating" ).as_boolean();
}

MirrorPlasma::MirrorPlasma( toml::value const& plasmaConfig )
	: pVacuumConfig( std::make_shared<VacuumMirrorConfiguration>( plasmaConfig ) )
{

	double TiTe = toml::find_or<double>( plasmaConfig, "IonToElectronTemperatureRatio", 1.0 );
	if ( TiTe <= 0.0 )
	{
		throw std::invalid_argument( "Ion to Electron temperature ratio must be a positive number" );
	}

	// Default to pure plasma
	Zeff = toml::find_or<double>( plasmaConfig, "Zeff", 1.0 );
	if ( Zeff <= 0.0 ) {
		throw std::invalid_argument( "Effective charge (Z_eff) must be positive!" );
	}

	try {
		ElectronDensity = toml::find<double>( plasmaConfig, "ElectronDensity" );
	} catch ( std::out_of_range &e ) {
		throw std::invalid_argument( "You must specify the electron density (ElectronDensity) in the [plasma] block" );
	}

	if ( ElectronDensity <= 0.0 ) {
		throw std::invalid_argument( "Electron density must be positive!" );
	}

	IonDensity = ElectronDensity / pVacuumConfig->IonSpecies.Charge;

	// If the temperature is -1.0, that indicates we are in a Temperature Solve run
	// and the temperature will be solved for.
	if ( plasmaConfig.count( "ElectronTemperature" ) == 1 ) {
		ElectronTemperature = toml::find_or<double>( plasmaConfig, "ElectronTemperature", -1.0 );
		IonTemperature = ElectronTemperature * TiTe;
	} else if ( plasmaConfig.count( "ElectronTemperature" ) == 0 ) {
		ElectronTemperature = -1.0;
		IonTemperature = -1.0;
	} else {
		throw std::invalid_argument( "Electron Temperature specified more than once!" );
	}

	if ( plasmaConfig.count( "NeutralDensity" ) == 1 ) {
		NeutralSource = 0;
		NeutralDensity = plasmaConfig.at( "NeutralDensity" ).as_floating();
	} else {
		NeutralSource = 0;
		NeutralDensity = 0;
	}

}

// From NRL Formulary p34
double MirrorPlasma::LogLambdaElectron() const
{
	return 17.;
	// Convert to NRL Formulary Units
	double neNRL = ElectronDensity * 1e14; // We use 10^20 / m^3 = 10^14 / cm^3
	double TeNRL = ElectronTemperature * 1000; // We normalize to 1 keV they use 1 ev
	// Convert log( ne^.5 Te^-1.25 ) = .5*log(ne) - 1.25 * log(Te) 
	double logTe = ::log( TeNRL );
	return 23.5 - 0.5*::log( neNRL ) + 1.25*logTe - ::sqrt( 1e-5 + ( logTe - 2 ) * ( logTe - 2 ) / 16 );
}

// tau_ee 
double MirrorPlasma::ElectronCollisionTime() const
{
	double tau_ee = ::sqrt( ElectronMass * ElectronTemperature * ReferenceTemperature / 2.0 ) * ( ElectronTemperature * ReferenceTemperature ) / 
	          ( ElectronDensity * ReferenceDensity * ::pow( GaussianElectronCharge, 4 ) * LogLambdaElectron() );
	return tau_ee;
}

// NRL Formulary p34
double MirrorPlasma::LogLambdaIon() const
{
	return 18.0;
	double TiNRL = IonTemperature * 1000;
	double niNRL = IonDensity * 1e14;
	double ZIon = pVacuumConfig->IonSpecies.Charge;
	return 23.0 - 0.5 * ::log( niNRL ) - 1.5 * ::log( ZIon * ZIon / TiNRL );
}

double MirrorPlasma::IonCollisionTime() const
{
	double PiThreeHalves = ::pow( M_PI, 1.5 ); // pi^(3/2)
	double TiThreeHalves = ::pow( IonTemperature * ReferenceTemperature, 1.5 );
	double ZIon = pVacuumConfig->IonSpecies.Charge;
	return 12 * ::sqrt( pVacuumConfig->IonSpecies.Mass * ProtonMass ) * PiThreeHalves * TiThreeHalves * VacuumPermittivity * VacuumPermittivity / ( ElectronDensity * ReferenceDensity * ::pow( ZIon * ElectronCharge, 4 ) * LogLambdaIon() );
}

// Leading order contribution to Phi_0 in O(M^2)
double MirrorPlasma::CentrifugalPotential() const
{
	return 0.5 * ( 1.0 - 1.0 / pVacuumConfig->MirrorRatio ) * MachNumber * MachNumber / ( pVacuumConfig->IonSpecies.Charge + 1 ); // IonTemperature/ElectronTemperature );
}

/* OLD CODE was in s^-1 not W!! */
double MirrorPlasma::ParallelElectronHeatLoss() const
{
	double R = pVacuumConfig->MirrorRatio;
	double tau_ee = ElectronCollisionTime();
	double x0 = CentrifugalPotential();
	// The integral in Pastukhov is, with x = T/e Phi
	// evaluates to 1 + sqrt(pi/x) * exp(x) * erfc( sqrt(x) )
	double PastukhovIntegral = 1.0 + ::sqrt( M_PI / x0 ) * ::exp( x0 ) * std::erfc( ::sqrt( x0 ) );
	double ElectronStoredEnergy = 1.5 * ElectronDensity * ElectronTemperature * ReferenceDensity * ReferenceTemperature;
	double LossRate = ( M_2_SQRTPI / tau_ee ) * ( 2.0 * R / ( 2.0 * R + 1 ) ) * ( 1.0 / ::log( 4.0 * R + 2.0 ) ) * ( 2./3. + ( 1./x0 ) * PastukhovIntegral ) * ::exp( - x0 );
	return LossRate * ElectronStoredEnergy;
}


double MirrorPlasma::ParallelKineticEnergyLoss() const
{
	double R = pVacuumConfig->MirrorRatio;
	double tau_ii = IonCollisionTime();
	double ZIon = pVacuumConfig->IonSpecies.Charge;
	double TiTe = IonTemperature / ElectronTemperature;

	// Centrifugal potential is normalised to e/T_e
	double x0 = CentrifugalPotential() * ZIon / TiTe;

	double PastukhovIntegral = 1 + ::sqrt( M_PI / x0 ) * ::exp( x0 ) * std::erfc( ::sqrt( x0 ) );

	// This is a bad estimate, as some deconfining will occur to ensure quasineutrality.
	double ParticleLossRate = ( M_2_SQRTPI / tau_ii ) * ( 2.0 * R / ( 2.0 * R + 1 ) ) * ( 1.0 / ::log( 4.0 * R + 2.0 ) ) * ( 2./3. + ( 1./x0 ) * PastukhovIntegral ) * ::exp( - x0 );

	double EndExpansionRatio = 1.0; // If particles are hitting the wall at a radius larger than they are confined at, this increases momentum loss per particle 
	double KineticEnergyPerParticle = 0.5 * ElectronTemperature * ReferenceTemperature * ( EndExpansionRatio * EndExpansionRatio ) *  MachNumber * MachNumber;
	// Convert to W/m^3 for consistency across loss rates
	return ParticleLossRate * ElectronDensity * ReferenceDensity * KineticEnergyPerParticle;
};


double MirrorPlasma::ClassicalIonHeatLoss() const
{
	double omega_ci = IonCyclotronFrequency();
	double IonMass = pVacuumConfig->IonSpecies.Mass * ProtonMass;
	double kappa_perp = 2.0 * IonTemperature * ReferenceTemperature / ( IonMass * omega_ci * omega_ci * IonCollisionTime() );
	double L_T = ( pVacuumConfig->PlasmaColumnWidth / 2.0 );
	double IonThermalEnergy = 1.5 * IonDensity * IonTemperature * ReferenceTemperature * ReferenceDensity;
	return kappa_perp * IonThermalEnergy /( L_T * L_T );
}

double MirrorPlasma::ClassicalHeatLosses() const 
{
	return ClassicalIonHeatLoss();
}

// Formula from the NRL formulary page 58 and the definition of Z_eff:
// n_e Z_eff = Sum_i n_i Z_i^2 (sum over all species that aren't electrons)
double MirrorPlasma::BremsstrahlungLosses() const
{
	// NRL formulary with reference values factored out
	// TODO: BUG HERE FOR COMPAT
	// Rteurn units are W/m^3
	return 169 * ::sqrt( 1000 * ElectronTemperature ) * Zeff * Zeff* ElectronDensity * ElectronDensity;
}

double MirrorPlasma::Beta() const
{
	// From Plasma Formulary, so convert to /cm^3 , eV, and Gauss
	double ne_Formulary = ElectronDensity * 1e-6; //TODO: BUG XXX
	double Te_Formulary = ElectronTemperature * 1e3;
	double ni_Formulary = IonDensity * ReferenceDensity * 1e-6;
	double Ti_Formulary = IonTemperature * 1e3;
	double MagField_Formulary = pVacuumConfig->CentralCellFieldStrength * 1e4;
	return 4.03e-11 * ( ne_Formulary * Te_Formulary + ni_Formulary * Ti_Formulary ) / ( MagField_Formulary * MagField_Formulary );
}

double MirrorPlasma::DebyeLength() const
{
	// l_debye^2 = epsilon_0 * (k_B T_e) / n_e * e^2
	double LambdaDebyeSquared = VacuumPermittivity * ( ElectronTemperature * ReferenceTemperature ) / ( ElectronDensity * ReferenceDensity * ElectronCharge * ElectronCharge );
	return ::sqrt( LambdaDebyeSquared );
}


double MirrorPlasma::NuStar() const
{
	return pVacuumConfig->PlasmaLength / ( IonCollisionTime() * SoundSpeed() );
}


double MirrorPlasma::KineticEnergy() const
{
	double IonMass = ProtonMass * pVacuumConfig->IonSpecies.Mass;
	double u = MachNumber * SoundSpeed();
	return .5 * IonMass * IonDensity * ReferenceDensity * u * u * pVacuumConfig->PlasmaVolume();
}

double MirrorPlasma::ThermalEnergy() const
{
	return 1.5 * ( ElectronDensity * ElectronTemperature + IonDensity * IonTemperature ) * ReferenceDensity * ReferenceTemperature * pVacuumConfig->PlasmaVolume(); 
}

// Braginskii eta_1
double MirrorPlasma::ClassicalViscosity() const
{
	double omega_ci = IonCyclotronFrequency();
	return ( 3.0 / 10.0 ) * ( pVacuumConfig->IonSpecies.Charge * ElectronDensity * ReferenceDensity * IonTemperature * ReferenceTemperature ) / ( omega_ci * omega_ci * IonCollisionTime() );
}

// Viscous heating = eta * u^2 / L_u^2
double MirrorPlasma::ViscousHeating() const
{
	double L_u = ( pVacuumConfig->PlasmaColumnWidth / 2.0 );
	double Velocity = MachNumber * SoundSpeed();
	double VelocityShear = Velocity / L_u;

	return ClassicalViscosity() * VelocityShear * VelocityShear;
}

double MirrorPlasma::AlfvenMachNumber() const
{
	// v_A^2 = B^2 / mu_0 * n_i * m_i
	// c_s^2 = T_e / m_i
	double MassDensity = IonDensity * ReferenceDensity * pVacuumConfig->IonSpecies.Mass * ProtonMass;
	double AlfvenSpeed = pVacuumConfig->CentralCellFieldStrength / ::sqrt( PermeabilityOfFreeSpace * MassDensity );

	return MachNumber * ( SoundSpeed() / AlfvenSpeed );
}

double MirrorPlasma::CollisionalTemperatureEquilibrationTime() const
{
	return IonCollisionTime()/( (3./pVacuumConfig->IonSpecies.Mass)*(ElectronMass/ProtonMass) );
}

double MirrorPlasma::IonToElectronHeatTransfer() const
{
	return 0.0;
}

double MirrorPlasma::IonHeatLosses() const
{
	return 2*ClassicalIonHeatLoss();
}

double MirrorPlasma::ElectronHeatLosses() const
{
	return 2*ParallelElectronHeatLoss() + BremsstrahlungLosses();
}

double MirrorPlasma::IonHeating() const
{
	double Heating = ViscousHeating();

	return Heating - IonToElectronHeatTransfer();
}

double MirrorPlasma::ElectronHeating() const
{
	double Heating = pVacuumConfig->AuxiliaryHeating;
	if ( pVacuumConfig->AlphaHeating ) {
		/*
		 * If the slowing-down time is longer than the energy confinement time
		 * reduce the alpha heating by that amount
		 */
		double AlphaHeatingFraction = EnergyConfinementTime() / SlowingDownTime();
		AlphaHeatingFraction = 1.0;
		// 1e6 as we are using W/m^3 and the formulary was in MW/m^3
		Heating += AlphaHeating() * 1e6 * AlphaHeatingFraction;

	}
	return Heating + IonToElectronHeatTransfer();
}


void MirrorPlasma::SetMachFromVoltage()
{
	// u = E x B / B^2 
	// M = u/c_s ~ (V/aB)/cs
	MachNumber = pVacuumConfig->ImposedVoltage / ( pVacuumConfig->PlasmaColumnWidth * pVacuumConfig->CentralCellFieldStrength * SoundSpeed() );
}


