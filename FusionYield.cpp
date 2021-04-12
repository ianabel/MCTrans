
#include "FusionYield.hpp"
#include <cmath>

// Energy density of alphas from fusion reactions in MW/m^3
double FusionAlphaPowerDensity( Plasma const& plasma, Configuration const& conf )
{
	double Ti = plasma.TiTe * plasma.ElectronTemperature;
	// From plasma formulary
	// Note n_e is in units of 10^20/m^3, Formulary is in particles per cc
	// we keep intermediate result in cm^3/s
	double SigmaVelocityAverage = 3.68e-12 * ::pow( Ti, -2./3. ) * ::exp( -19.94 * ::pow( Ti, -1./3. ) );
	// we assume n_D = n_T = .5 * n_e (could include impurity dilution?)
	double PowerDensity = 5.6e-13 * ( 0.5 * 0.5 * plasma.ElectronDensity * plasma.ElectronDensity * 1e28 ) * SigmaVelocityAverage;
	// Units are W/cc from formulary, which is conveniently also MW/m^3
	// Formulary only includes charged particle energy, so we include the neutron
	return PowerDensity;
}

// Neutron output in MW
double NeutronOutput( Plasma const& plasma, Configuration const& conf )
{
	return PlasmaVolume( conf ) * FusionAlphaPowerDensity( plasma, conf ) * ( 14.1 / 3.52 );
}

// Flux of neutrons through wall in MW/m^2
double NeutronWallLoading( Plasma const& plasma, Configuration const& conf )
{
	return NeutronOutput( plasma, conf )/WallArea( conf );
}

// MW Th 
double ThermalPowerOutput( Plasma const& plasma, Configuration const& conf )
{
	// Thermal power has several components:
	// 1) Heat leaving the plasma, which in steady state is Aux + Alphas
	double ExhaustHeat = conf.AuxiliaryHeating + FusionAlphaPowerDensity( plasma, conf ) * PlasmaVolume( conf );
	// 2)  Neutrons absorbed in the blanket
	double NeutronYield = NeutronOutput( plasma, conf );
	// 3) Heat Generated in the blanket from breeding reactions (which is 4.8MeV per triton, which is 4.8MeV per incident neutron assuming a TBR of 1), and thus
	double BreedingYield = NeutronYield * ( 4.8 / 14.1 );
	// Note that this is not a negligible amount of heat.
	return ExhaustHeat + NeutronYield + BreedingYield;
}

// in neutrons / s
double DDNeutronRate( Plasma const& plasma, Configuration const& conf )
{
	double Ti = plasma.TiTe * plasma.ElectronTemperature;
	// From plasma formulary
	double SigmaVelocityAverage = 2.33e-14 * ::pow( Ti, -2./3. ) * ::exp( -18.76 * ::pow( Ti, -1./3. ) );

	// n_e is in units of 10^20/m^3, Formulary is in particles per cc
	// we assume n_D = n_e
	double nD_cc = plasma.ElectronDensity * 1e14;
	double NeutronRatePerCC = nD_cc * nD_cc * SigmaVelocityAverage;

	return NeutronRatePerCC * PlasmaVolume( conf ) * 1e6;
}



