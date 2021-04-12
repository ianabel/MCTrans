#include <iostream>
#include <iomanip>
#include <cmath>

#include "Centrifugal.hpp"
#include "Config.hpp"
#include "FusionYield.hpp"
#include "AlphaHeating.hpp"

void PrintWithUnit( double value, std::string const& unit )
{
	if ( ::fabs( value ) > 1e9 || ::fabs( value ) < 1e3 )
	{
		std::cout << value << " " << unit;
		return;
	}

	if ( ::fabs( value ) > 1e3 && ::fabs( value ) <= 1e6 )
		std::cout << value/1e3 << " k" << unit;
	else if ( ::fabs( value ) > 1e6 && ::fabs( value ) <= 1e9 )
		std::cout << value/1e6 << " M" << unit;

	return;
}
	


void PrintReport( Plasma const& plasma, Configuration &conf )
{
	std::cout.precision( 3 );



	std::cout << "The plasma is made up of electrons, " << plasma.Fuel;
	if ( plasma.LumpedImpurity )
		std::cout << " and lumped impurities" << std::endl;
	else if ( plasma.Zeff != 0.0 )
		std::cout << " and trace, radiating, impurities" << std::endl;
	else 
		std::cout << std::endl;

	std::cout << "Electron Density is " << plasma.ElectronDensity*ReferenceDensity << " m^-3" << std::endl;
	std::cout << "Ion Density is " << plasma.IonDensity * ReferenceDensity << " m^-3" << std::endl;
	if ( plasma.LumpedImpurity )
		std::cout << "Impurity Density is " << plasma.ImpurityDensity * ReferenceDensity << " m^-3" << std::endl;
	std::cout << std::endl;
	std::cout << "Electron Temperature is " << plasma.ElectronTemperature << " keV" << std::endl;
	std::cout << "Ion Temperature is " << plasma.TiTe * plasma.ElectronTemperature << " keV" << std::endl;

	std::cout << "Radius of central cell (first wall) is " << conf.WallRadius << " m" << std::endl;
	std::cout << "\t Inner radius of the plasma is " << conf.AxialGapDistance << " m" << std::endl;
	std::cout << "\t Outer radius of the plasma is " << conf.PlasmaColumnWidth + conf.AxialGapDistance << " m" << std::endl;

	std::cout << "Outer Radius of the plasma at the mirror throat is " << (conf.PlasmaColumnWidth + conf.AxialGapDistance) /::sqrt( conf.MirrorRatio ) << " m" << std::endl;
	std::cout << std::endl;
	std::cout << "Axial Length of the plasma is " << conf.PlasmaLength << " m" << std::endl;
	std::cout << "Plasma Volume is " << PlasmaVolume( conf ) << " m^3" << std::endl;
	std::cout << "" << std::endl;
	std::cout << "Magnetic Field in the Central Cell is " << conf.CentralCellFieldStrength << " T" << std::endl;
	std::cout << "Magnetic Field at the Mirror Throat is " << conf.MirrorRatio * conf.CentralCellFieldStrength << " T" << std::endl;
	std::cout << "" << std::endl;
	std::cout << "Ion Larmor Radius in the central cell is " << IonLarmorRadius( plasma, conf ) << " m " << std::endl;
	std::cout << "Typical plasma scale lenghts are ~ " << conf.PlasmaColumnWidth/2.0 << " m  = " << conf.PlasmaColumnWidth / ( 2.0 * IonLarmorRadius( plasma, conf ) ) << " rho_i " << std::endl;

	std::cout << std::endl;

	if ( conf.AuxiliaryHeating > 0 ) {
		std::cout << "Auxiliary Heating of " << conf.AuxiliaryHeating * 1e3  << " kW was included" << std::endl;
		std::cout << "" << std::endl;
	} else {
		std::cout << "No auxiliary heating was included in this calculation." << std::endl;
	}

	if ( conf.IncludeAlphaHeating ) {
		std::cout << " Self-consistent Alpha Heating was included in this calculation. " << std::endl;

		std::cout << " Alphas provide " << AlphaHeating( plasma, conf )*PlasmaVolume( conf ) << " MW of heating" << std::endl;
		std::cout << " Prompt Alpha Losses towards the end plates give " << AlphaPromptLosses( plasma, conf ) * PlasmaVolume( conf ) << " MW of energy losses" << std::endl;
		double TauSD = SlowingDownTime( plasma );
		std::cout << " The alpha particle slowing-down time is " << TauSD << " s" << std::endl;
	}

	std::cout << std::endl;
	std::cout << "Operating Mach number is " << conf.MachNumber << std::endl;
	std::cout << "Alfven Mach number is " << AlfvenMachNumber( plasma, conf ) << std::endl;
	std::cout << std::endl;

	std::cout << "Viscous Heating is ";
	PrintWithUnit( ViscousHeating( plasma, conf ) * PlasmaVolume( conf ),"W" ); std::cout << std::endl;
	if ( plasma.Zeff > 0.0 )
	{
		std::cout << "Bremsstrahlung losses are ";
		PrintWithUnit( BremsstrahlungLosses( plasma )*PlasmaVolume( conf ), "W" ); std::cout << std::endl;
	}
	else
		std::cout << "No radiation losses were included." << std::endl;
	std::cout << std::endl;
		
	std::cout << "Total potential drop is ";
	PrintWithUnit( ElectricPotential( plasma, conf ), "V" );
	std::cout << std::endl;

	double ViscousHeatingRate = ViscousHeating( plasma, conf ) * PlasmaVolume( conf );
	double ParallelMomentumLoss = ParallelKineticEnergyLoss( plasma, conf ) * PlasmaVolume( conf );
	std::cout << "Power Required (at the plasma) to support rotation ";
	PrintWithUnit( ViscousHeatingRate + ParallelMomentumLoss, "W" );
	std::cout << std::endl;
	std::cout << std::endl;

	double KineticStoredEnergy = KineticEnergy( plasma, conf );
	double ThermalStoredEnergy = ThermalEnergy( plasma, conf );
	std::cout << "Total Stored Energy is "; PrintWithUnit( KineticStoredEnergy + ThermalStoredEnergy, "J" ); std::cout << std::endl;
	std::cout << " of which "; PrintWithUnit( KineticStoredEnergy, "J" ); std::cout << " is kinetic energy in the rotation" << std::endl;
	std::cout << "      and "; PrintWithUnit( ThermalStoredEnergy, "J" ); std::cout << " is thermal energy of the plasma" << std::endl;


	std::cout << "Energy Confinement Time is " << 1.0/TauEInverse( plasma, conf ) << " s" << std::endl;
	std::cout << "Of which" << std::endl;
	std::cout << "\tConfinement time from parallel losses is " << 1.0/ParallelElectronHeatLoss( plasma, conf ) << " s" << std::endl;
	std::cout << "\tConfinement time from perpendicular losses is " << 1.0/ClassicalHeatLossRate( plasma, conf ) << " s" << std::endl;
	std::cout << std::endl;

	double IonElectronEquilibrationTime = IonCollisionTime( plasma )/( (3./plasma.Mu)*(ElectronMass/ProtonMass) );
	std::cout << "Ion-Electron Temperature Equilibration Time is " << IonElectronEquilibrationTime << " s" << std::endl;
	double EquilibrationRatio = IonElectronEquilibrationTime * TauEInverse( plasma, conf );
	/*
	if ( EquilibrationRatio > 0.9 )
		std::cout << "WARNING: Assumption of fixed temperature ratio may be invalid." << std::endl;
	*/


	std::cout << std::endl;
	std::cout << "Dimensionless parameters:" << std::endl;
	std::cout << "\t β  = " << Beta( plasma, conf ) * 100 << "%" << std::endl;
	std::cout << "\t ν* = " << NuStar( plasma, conf ) << " (ions) " << std::endl; 
	std::cout << "\t ρ* = " << ( 2.0 * IonLarmorRadius( plasma, conf ) ) / conf.PlasmaColumnWidth << std::endl; 
	

	std::cout << std::endl;
	std::cout << "Fusion Triple Product " << std::endl;
	double TripleProduct = plasma.IonDensity * ReferenceDensity * plasma.ElectronTemperature * plasma.TiTe / TauEInverse( plasma, conf );
	std::cout << "\t n T τ = " << TripleProduct << " keV s /m^3" << std::endl;
	
	std::cout << std::endl;



	if ( conf.ReportNuclear ) {
		std::cout << " === Nuclear Reactions Assuming Deuterium Fuel === " << std::endl;
		std::cout << " Neutrons per second: " << DDNeutronRate( plasma, conf ) << std::endl;
		std::cout << std::endl;
		std::cout << " === Reactor Output Assuming D/T Reactor Fuel === " << std::endl;

		double FusionAlphaPower = FusionAlphaPowerDensity( plasma, conf )*PlasmaVolume( conf );
		double FusionNeutronPower = ( 14.1/3.52 ) * FusionAlphaPower;
		std::cout << "Fusion Power Output as Alphas:   " << FusionAlphaPower   << " MW" << std::endl;
		std::cout << "                    as Neutrons: " << FusionNeutronPower << " MW" << std::endl;
		std::cout << "Neutron Wall Loading is " << NeutronWallLoading( plasma, conf ) << " MW/m^2" << std::endl;
		std::cout << "Total Thermal Power Output is " << ThermalPowerOutput( plasma, conf ) << " MW" << std::endl;
	}

}

int main( int argc, char** argv )
{
	Plasma plasma;
	Configuration conf;

	std::string fname( "Mirror.conf" );
	if ( argc == 2 )
		fname = argv[ 1 ];
	if ( argc > 2 )
	{
		std::cerr << "Usage: MCTrans++ ConfigFile.conf" << std::endl;
		return 1;
	}

	LoadConfig( fname, plasma, conf );

	if ( conf.MachSolve )
		SetMachNumber( plasma, conf );
	else if ( conf.TempSolve )
		SetTemperature( plasma, conf );

	PrintReport( plasma, conf );
	return 0.0;

}
