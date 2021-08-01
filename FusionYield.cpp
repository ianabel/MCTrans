
#include "FusionYield.hpp"
#include "MirrorPlasma.hpp"
#include <cmath>

double FusionReactions::SigmaAverageDT( double Ti ) {
	// From plasma formulary
	// Note n_e is in units of 10^20/m^3, Formulary is in particles per cc
	// we keep intermediate result in cm^3/s
	return 3.68e-12 * ::pow( Ti, -2./3. ) * ::exp( -19.94 * ::pow( Ti, -1./3. ) );
}

double FusionReactions::SigmaAverageDD( double Ti ) {
	// From plasma formulary
	// Note n_e is in units of 10^20/m^3, Formulary is in particles per cc
	// we keep intermediate result in cm^3/s
	return 2.33e-14 * ::pow( Ti, -2./3. ) * ::exp( -18.76 * ::pow( Ti, -1./3. ) );
}

// Energy density of alphas from fusion reactions in MW/m^3
double MirrorPlasma::FusionAlphaPowerDensity() const
{
	// From plasma formulary
	// Note n_e is in units of 10^20/m^3, Formulary is in particles per cc
	// we keep intermediate result in cm^3/s
	double SigmaVelocityAverage = FusionReactions::SigmaAverageDT( IonTemperature );
	// we assume n_D = n_T = .5 * n_e (could include impurity dilution?)
	double nD_cc = .5 * IonDensity * 1e14;
	double nT_cc = .5 * IonDensity * 1e14;
	double PowerDensity = 5.6e-13 * nD_cc * nT_cc * SigmaVelocityAverage;
	// Units are W/cc from formulary, which is conveniently also MW/m^3
	return PowerDensity;
}

// Neutron output in MW
double MirrorPlasma::NeutronOutput() const
{
	return pVacuumConfig->PlasmaVolume() * FusionAlphaPowerDensity() * ( 14.1 / 3.52 );
}

// Flux of neutrons through wall in MW/m^2
double MirrorPlasma::NeutronWallLoading() const
{
	return NeutronOutput()/pVacuumConfig->WallArea();
}

// MW Th 
double MirrorPlasma::ThermalPowerOutput() const
{
	// Thermal power has several components:
	// 1) Heat leaving the plasma, which in steady state is Aux + Alphas
	double ExhaustHeat = pVacuumConfig->AuxiliaryHeating + FusionAlphaPowerDensity() * pVacuumConfig->PlasmaVolume();
	// 2)  Neutrons absorbed in the blanket
	double NeutronYield = NeutronOutput();
	// 3) Heat Generated in the blanket from breeding reactions (which is 4.8MeV per triton, which is 4.8MeV per incident neutron assuming a TBR of 1), and thus
	double BreedingYield = NeutronYield * ( 4.8 / 14.1 );
	// Note that this is not a negligible amount of heat.
	return ExhaustHeat + NeutronYield + BreedingYield;
}

// in neutrons / s
double MirrorPlasma::DDNeutronRate() const
{
	// From plasma formulary
	double SigmaVelocityAverage = FusionReactions::SigmaAverageDD( IonTemperature );

	// n_e is in units of 10^20/m^3, Formulary is in particles per cc
	// we assume n_D = n_e
	double nD_cc = IonDensity * 1e14;
	double NeutronRatePerCC = nD_cc * nD_cc * SigmaVelocityAverage;

	return NeutronRatePerCC * pVacuumConfig->PlasmaVolume() * 1e6;
}


double MirrorPlasma::AlphaProductionRate() const
{
	double SigmaVelocityAverage = FusionReactions::SigmaAverageDT( IonTemperature );
	double nD_cc = .5 * IonDensity * 1e14;
	double nT_cc = .5 * IonDensity * 1e14;
	// Returns alphas / m^3
	return nD_cc * nT_cc * SigmaVelocityAverage * 1e6;
}

