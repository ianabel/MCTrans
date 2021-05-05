
#include "MirrorPlasma.hpp"
#include <iostream>

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

void MirrorPlasma::PrintReport() const
{
	std::cout.precision( 3 );

	std::cout << "The plasma is made up of electrons, " << pVacuumConfig->IonSpecies.Name; 
	if ( Zeff > 0.0 )
		std::cout << " and trace, radiating, impurities" << std::endl;
	else 
		std::cout << std::endl;

	std::cout << "Electron Density is " << ElectronDensity*ReferenceDensity << " m^-3" << std::endl;
	std::cout << "Ion Density is " << IonDensity * ReferenceDensity << " m^-3" << std::endl;
	std::cout << std::endl;
	std::cout << "Electron Temperature is "; PrintWithUnit( ElectronTemperature * 1e3,"eV"); std::cout << std::endl;
	std::cout << "Ion Temperature is "; PrintWithUnit( IonTemperature * 1e3,"eV"); std::cout << std::endl;

	std::cout << "Radius of central cell (first wall) is " << pVacuumConfig->WallRadius << " m" << std::endl;

	std::cout << "\t Inner radius of the plasma is " << pVacuumConfig->AxialGapDistance << " m" << std::endl;
	std::cout << "\t Outer radius of the plasma is " << pVacuumConfig->PlasmaColumnWidth + pVacuumConfig->AxialGapDistance << " m" << std::endl;

	std::cout << "Outer Radius of the plasma at the mirror throat is " << (pVacuumConfig->PlasmaColumnWidth + pVacuumConfig->AxialGapDistance) /::sqrt( pVacuumConfig->MirrorRatio ) << " m" << std::endl;
	std::cout << std::endl;
	std::cout << "Axial Length of the plasma is " << pVacuumConfig->PlasmaLength << " m" << std::endl;
	std::cout << "Plasma Volume is " << pVacuumConfig->PlasmaVolume() << " m^3" << std::endl;
	std::cout << "" << std::endl;
	std::cout << "Magnetic Field in the Central Cell is " << pVacuumConfig->CentralCellFieldStrength << " T" << std::endl;
	std::cout << "Magnetic Field at the Mirror Throat is " << pVacuumConfig->MirrorRatio * pVacuumConfig->CentralCellFieldStrength << " T" << std::endl;
	std::cout << "" << std::endl;
	std::cout << "Ion Larmor Radius in the central cell is " << IonLarmorRadius() << " m " << std::endl;
	std::cout << "Typical plasma scale lenghts are ~ " << pVacuumConfig->PlasmaColumnWidth/2.0 << " m  = " << pVacuumConfig->PlasmaColumnWidth / ( 2.0 * IonLarmorRadius() ) << " rho_i " << std::endl;

	std::cout << std::endl;

	if ( pVacuumConfig->AuxiliaryHeating > 0 ) {
		std::cout << "Auxiliary Heating of " << pVacuumConfig->AuxiliaryHeating * 1e3  << " kW was included" << std::endl;
		std::cout << "" << std::endl;
	} else {
		std::cout << "No auxiliary heating was included in this calculation." << std::endl;
	}

	if ( pVacuumConfig->AlphaHeating ) {
		std::cout << " Self-consistent Alpha Heating was included in this calculation. " << std::endl;

		std::cout << " Alphas provide " << AlphaHeating()*pVacuumConfig->PlasmaVolume() << " MW of heating" << std::endl;
		std::cout << " Prompt Alpha Losses towards the end plates give " << AlphaPromptLosses() * pVacuumConfig->PlasmaVolume() << " MW of energy losses" << std::endl;
		double TauSD = SlowingDownTime();
		std::cout << " The alpha particle slowing-down time is " << TauSD << " s" << std::endl;
	}

	std::cout << std::endl;
	std::cout << "Operating Mach number is " << MachNumber << std::endl;
	std::cout << "Alfven Mach number is " << AlfvenMachNumber() << std::endl;
	std::cout << std::endl;

	std::cout << "Viscous Heating is ";
	PrintWithUnit( ViscousHeating() * pVacuumConfig->PlasmaVolume(),"W" ); std::cout << std::endl;
	if ( Zeff > 0.0 )
	{
		std::cout << "Bremsstrahlung losses are ";
		PrintWithUnit( BremsstrahlungLosses()*pVacuumConfig->PlasmaVolume(), "W" ); std::cout << std::endl;
	}
	else
		std::cout << "No radiation losses were included." << std::endl;
	std::cout << std::endl;
		
	std::cout << "Total potential drop is ";
	PrintWithUnit( ElectricPotential(), "V" );
	std::cout << std::endl;

	double ViscousHeatingRate = ViscousHeating() * pVacuumConfig->PlasmaVolume();
	double ParallelMomentumLoss = ParallelKineticEnergyLoss() * pVacuumConfig->PlasmaVolume();
	std::cout << "Power Required (at the plasma) to support rotation ";
	PrintWithUnit( ViscousHeatingRate + ParallelMomentumLoss, "W" );
	std::cout << std::endl;
	std::cout << std::endl;

	double KineticStoredEnergy = KineticEnergy();
	double ThermalStoredEnergy = ThermalEnergy();
	std::cout << "Total Stored Energy is "; PrintWithUnit( KineticStoredEnergy + ThermalStoredEnergy, "J" ); std::cout << std::endl;
	std::cout << " of which "; PrintWithUnit( KineticStoredEnergy, "J" ); std::cout << " is kinetic energy in the rotation" << std::endl;
	std::cout << "      and "; PrintWithUnit( ThermalStoredEnergy, "J" ); std::cout << " is thermal energy of the plasma" << std::endl;


	std::cout << "Energy Confinement Time is " << EnergyConfinementTime() << " s" << std::endl;
	std::cout << "Of which" << std::endl;
	
	std::cout << "\tConfinement time from parallel losses is " << ( 1.5  * ElectronDensity * ElectronTemperature * ReferenceDensity * ReferenceTemperature )/ParallelElectronHeatLoss() << " s" << std::endl;
	std::cout << "\tConfinement time from perpendicular losses is " << ( 1.5  * IonDensity * IonTemperature * ReferenceDensity * ReferenceTemperature )/ClassicalIonHeatLoss() << " s" << std::endl;
	std::cout << std::endl;

	double IonElectronEquilibrationTime = CollisionalTemperatureEquilibrationTime();
	std::cout << "Ion-Electron Temperature Equilibration Time is " << IonElectronEquilibrationTime << " s" << std::endl;
	double EquilibrationRatio = IonElectronEquilibrationTime / EnergyConfinementTime(); 
	/*
	if ( EquilibrationRatio > 0.9 )
		std::cout << "WARNING: Assumption of fixed temperature ratio may be invalid." << std::endl;
	*/


	std::cout << std::endl;
	std::cout << "Dimensionless parameters:" << std::endl;
	std::cout << "\t β  = " << Beta() * 100 << "%" << std::endl;
	std::cout << "\t ν* = " << NuStar() << " (ions) " << std::endl; 
	std::cout << "\t ρ* = " << ( 2.0 * IonLarmorRadius() ) / pVacuumConfig->PlasmaColumnWidth << std::endl; 
	

	std::cout << std::endl;
	std::cout << "Fusion Triple Product " << std::endl;
	double TripleProduct = IonDensity * ReferenceDensity * IonTemperature * EnergyConfinementTime();
	std::cout << "\t n T τ = " << TripleProduct << " keV s /m^3" << std::endl;
	
	std::cout << std::endl;



	if ( pVacuumConfig->ReportNuclearDiagnostics ) {
		std::cout << " === Nuclear Reactions Assuming Deuterium Fuel === " << std::endl;
		std::cout << " Neutrons per second: " << DDNeutronRate() << std::endl;
		std::cout << std::endl;
		std::cout << " === Reactor Output Assuming D/T Reactor Fuel === " << std::endl;

		double FusionAlphaPower = FusionAlphaPowerDensity()*pVacuumConfig->PlasmaVolume();
		double FusionNeutronPower = ( 14.1/3.52 ) * FusionAlphaPower;
		std::cout << "Fusion Power Output as Alphas:   " << FusionAlphaPower   << " MW" << std::endl;
		std::cout << "                    as Neutrons: " << FusionNeutronPower << " MW" << std::endl;
		std::cout << "Neutron Wall Loading is " << NeutronWallLoading() << " MW/m^2" << std::endl;
		std::cout << "Total Thermal Power Output is " << ThermalPowerOutput() << " MW" << std::endl;
	}

}
