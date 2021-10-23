
#include "MirrorPlasma.hpp"
#include <iostream>

void PrintWithUnit( std::ostream& out, double value, std::string const& unit )
{

	if ( ::fabs( value ) > 1e3 && ::fabs( value ) <= 1e6 )
		out << value/1e3 << " k" << unit;
	else if ( ::fabs( value ) > 1e6 && ::fabs( value ) <= 1e9 )
		out << value/1e6 << " M" << unit;
	else if ( ::fabs( value ) < 1 && ::fabs( value ) >= 1e-3 )
		out << value*1e3 << " m" << unit;
	else if ( ::fabs( value ) < 1e-3 && ::fabs( value ) >= 1e-6 )
		out << value*1e6 << " µ" << unit;
	else if ( ::fabs( value ) < 1e-6 && ::fabs( value ) >= 1e-9 )
		out << value*1e9 << " p" << unit;
	else
		out << value << " " << unit;

	return;
}


void MirrorPlasma::PrintReport()
{

	std::ostream *p_out;
	if ( pVacuumConfig->OutputFile == "" )
		p_out = &std::cout;
	else
	{
		p_out = new std::fstream( pVacuumConfig->OutputFile, std::ios_base::out | std::ios_base::trunc );
	}

	std::ostream &out = *p_out;

	out.precision( 3 );

	out << "The plasma is made up of electrons, " << pVacuumConfig->IonSpecies.Name; 
	if ( Zeff > 0.0 )
		out << " and trace, radiating, impurities" << std::endl;
	else 
		out << std::endl;

	out << "Electron Density is " << ElectronDensity*ReferenceDensity << " m^-3" << std::endl;
	out << "Ion Density is " << IonDensity * ReferenceDensity << " m^-3" << std::endl;
	out << std::endl;
	out << "Electron Temperature is "; PrintWithUnit( out, ElectronTemperature * 1e3,"eV"); out << std::endl;
	out << "Ion Temperature is "; PrintWithUnit( out, IonTemperature * 1e3,"eV"); out << std::endl;

	out << "Radius of central cell (first wall) is " << pVacuumConfig->WallRadius << " m" << std::endl;

	out << "\t Inner radius of the plasma is " << pVacuumConfig->AxialGapDistance << " m" << std::endl;
	out << "\t Outer radius of the plasma is " << pVacuumConfig->PlasmaColumnWidth + pVacuumConfig->AxialGapDistance << " m" << std::endl;

	out << "Outer Radius of the plasma at the mirror throat is " << (pVacuumConfig->PlasmaColumnWidth + pVacuumConfig->AxialGapDistance) /::sqrt( pVacuumConfig->MirrorRatio ) << " m" << std::endl;
	out << std::endl;
	out << "Axial Length of the plasma is " << pVacuumConfig->PlasmaLength << " m" << std::endl;
	out << "Plasma Volume is " << pVacuumConfig->PlasmaVolume() << " m^3" << std::endl;
	out << "" << std::endl;
	out << "Magnetic Field in the Central Cell is " << pVacuumConfig->CentralCellFieldStrength << " T" << std::endl;
	out << "Magnetic Field at the Mirror Throat is " << pVacuumConfig->MirrorRatio * pVacuumConfig->CentralCellFieldStrength << " T" << std::endl;
	out << "" << std::endl;
	out << "Ion Larmor Radius in the central cell is " << IonLarmorRadius() << " m " << std::endl;
	out << "Typical plasma scale lenghts are ~ " << pVacuumConfig->PlasmaColumnWidth/2.0 << " m  = " << pVacuumConfig->PlasmaColumnWidth / ( 2.0 * IonLarmorRadius() ) << " rho_i " << std::endl;

	out << std::endl;

	if ( pVacuumConfig->AuxiliaryHeating > 0 ) {
		out << "Auxiliary Heating of " << pVacuumConfig->AuxiliaryHeating * 1e3  << " kW was included" << std::endl;
		out << "" << std::endl;
	} else {
		out << "No auxiliary heating was included in this calculation." << std::endl;
		out << std::endl;
	}

	if ( pVacuumConfig->AlphaHeating ) {
		out << "Self-consistent Alpha Heating was included in this calculation. " << std::endl;

		out << "  Alphas provide " << AlphaHeating()*pVacuumConfig->PlasmaVolume() << " MW of heating" << std::endl;
		out << "  Alpha particles are lost directly at a rate of " << PromptAlphaLossFraction() * AlphaProductionRate() * pVacuumConfig->PlasmaVolume() << " /s" << std::endl;
		out << "  Prompt Alpha Losses towards the end plates give " << AlphaPromptLosses() * pVacuumConfig->PlasmaVolume() << " MW of energy losses" << std::endl;

		double TauSD = SlowingDownTime();
		out << "  The alpha particle slowing-down time is " << TauSD << " s" << std::endl;
	}

	out << std::endl;
	out << "Operating Mach number is " << MachNumber << std::endl;
	out << "Alfven Mach number is " << AlfvenMachNumber() << std::endl;
	out << std::endl;
	double v = SoundSpeed() * MachNumber;
	out << "Velocity is " << v << " m/s" << std::endl;
	double Rmid = pVacuumConfig->PlasmaCentralRadius();
	out << "Angular Velocity at the plasma centre (R = " << Rmid << " m) is " << v/Rmid << " /s" << std::endl;
	out << std::endl;

	out << "Viscous Heating is ";
	PrintWithUnit( out, ViscousHeating() * pVacuumConfig->PlasmaVolume(),"W" ); out << std::endl;
	if ( Zeff > 0.0 )
	{
		out << "Bremsstrahlung losses are ";
		PrintWithUnit( out, BremsstrahlungLosses()*pVacuumConfig->PlasmaVolume(), "W" ); out << std::endl;
	}
	else
		out << "No radiation losses were included." << std::endl;
	out << std::endl;
		
	out << "Total potential drop is ";
	PrintWithUnit( out, ElectricPotential(), "V" );
	out << std::endl;

	double JRadial = ::fabs( RadialCurrent() );
	out << "Radial Current Drawn from Power Supply "; PrintWithUnit( out, JRadial, "A" ); out << std::endl;
	out << "Power Required (at the plasma) to support rotation ";
	PrintWithUnit( out, ElectricPotential() * JRadial, "W" );
	out << std::endl;
	double omega = v/Rmid;
	out << "\t Power Loss from viscous torque "; PrintWithUnit( out, ViscousTorque()*omega*pVacuumConfig->PlasmaVolume(), "W" ); out << std::endl;
	out << "\t Power Loss from parallel loss  "; PrintWithUnit( out, ParallelAngularMomentumLossRate()*omega*pVacuumConfig->PlasmaVolume(), "W" ); out << std::endl;
	out << std::endl;

	double KineticStoredEnergy = KineticEnergy();
	double ThermalStoredEnergy = ThermalEnergy();
	out << "Total Stored Energy is "; PrintWithUnit( out, KineticStoredEnergy + ThermalStoredEnergy, "J" ); out << std::endl;
	out << " of which "; PrintWithUnit( out, KineticStoredEnergy, "J" ); out << " is kinetic energy in the rotation" << std::endl;
	out << "      and "; PrintWithUnit( out, ThermalStoredEnergy, "J" ); out << " is thermal energy of the plasma" << std::endl;


	out << "Energy Confinement Time is "; PrintWithUnit( out, EnergyConfinementTime(),"s" ); out << std::endl;
	out << "Of which" << std::endl;
	
	double ParallelConfinementTime = ( 1.5  * ElectronDensity * ElectronTemperature + 1.5 * IonDensity * IonTemperature ) * (  ReferenceDensity * ReferenceTemperature )/( ParallelElectronHeatLoss() + ParallelIonHeatLoss() );
	double PerpConfinementTime = ( 1.5  * IonDensity * IonTemperature * ReferenceDensity * ReferenceTemperature )/ClassicalIonHeatLoss();
	out << "\tConfinement time from parallel losses is "; PrintWithUnit( out, ParallelConfinementTime, "s" ); out << std::endl;
	out << "\tConfinement time from perpendicular losses is "; PrintWithUnit( out, PerpConfinementTime, "s" ); out << std::endl;
	out << std::endl;

	double ParallelParticleConfinementTime = ( IonDensity * ReferenceDensity ) /(  ParallelIonParticleLoss() );
	out << "Particle (Ion) Confinement Time ~= "; PrintWithUnit( out, ParallelParticleConfinementTime, "s" ); out << std::endl;
	out << std::endl;

	out << "Ion-Electron Temperature Equilibration Time is "; PrintWithUnit( out, CollisionalTemperatureEquilibrationTime(),"s" ); out << std::endl;

	/*
	out << std::endl;
	out << "Neutral gas must be provided at a rate of " << NeutralSource * pVacuumConfig->PlasmaVolume() << " particles/s to refuel the plasma" << std::endl;
	out << "This leads to a steady-state neutral density of " << NeutralDensity * ReferenceDensity << "/m^3" << std::endl;
	*/

	out << std::endl;
	double ElectronLossRate = ParallelElectronParticleLoss() + ClassicalElectronParticleLosses();
	out << "Plasma must be provided at a rate of " << ElectronLossRate*pVacuumConfig->PlasmaVolume() << " electrons /s to maintain steady-state" << std::endl;

	out << std::endl;
	double Resistance = ElectricPotential() / JRadial;
	double Capacitance = 2.0*KineticStoredEnergy / ( ElectricPotential() * ElectricPotential() );
	out << "Electrical Properties of the Plasma:" << std::endl;
	out << "\tResistance  = "; PrintWithUnit( out, Resistance,  "Ω" ); out << std::endl;
	out << "\tCapacitance = "; PrintWithUnit( out, Capacitance, "F" ); out << std::endl;



	out << std::endl;
	out << "Dimensionless parameters:" << std::endl;
	out << "\t β  = " << Beta() * 100 << "%" << std::endl;
	out << "\t ν* = " << NuStar() << " (ions) " << std::endl; 
	out << "\t ρ* = " << ( 2.0 * IonLarmorRadius() ) / pVacuumConfig->PlasmaColumnWidth << std::endl; 
	out << "\t Omega_i tau_ii = " << IonCyclotronFrequency()*IonCollisionTime() << std::endl; 

	

	out << std::endl;
	out << "Fusion Triple Product " << std::endl;
	double TripleProduct = IonDensity * ReferenceDensity * IonTemperature * EnergyConfinementTime();
	out << "\t n T τ = " << TripleProduct << " keV s /m^3" << std::endl;
	
	out << std::endl;

	if ( pVacuumConfig->ReportThrust ) {
			out << "Ions lost from the central cell provide " << ParallelIonThrust() << " N of thrust" << std::endl;
			if ( pVacuumConfig->AlphaHeating )
				out << "Prompt Alpha losses contribute " << PromptAlphaThrust() << " N of thrust" << std::endl;
	}


	if ( pVacuumConfig->ReportNuclearDiagnostics ) {
		if ( pVacuumConfig->IonSpecies.Name == "Deuterium" ) {
			out << " === Nuclear Reactions Assuming Deuterium Fuel === " << std::endl;
			out << " D/D Neutrons (2.45 MeV) per second: " << DDNeutronRate() << std::endl;
			// out << " D/T Neutrons (14.1 MeV) from Fusion-Produced-Tritium: " << DDEnergeticNeutronRate() << std::endl;
			out << std::endl;
		} else if ( pVacuumConfig->IonSpecies.Name == "Deuterium/Tritium Fuel" ) {
			out << " === Reactor Output Assuming D/T Reactor Fuel === " << std::endl;

			double FusionAlphaPower = FusionAlphaPowerDensity()*pVacuumConfig->PlasmaVolume();
			double FusionNeutronPower = ( 14.1/3.52 ) * FusionAlphaPower;
			out << "Fusion Power Output as Alphas:   " << FusionAlphaPower   << " MW" << std::endl;
			out << "                    as Neutrons: " << FusionNeutronPower << " MW" << std::endl;
			out << "Neutron Wall Loading is " << NeutronWallLoading() << " MW/m^2" << std::endl;
			out << "Total Thermal Power Output is " << ThermalPowerOutput() << " MW" << std::endl;

		} else {
			out << "No Nuclear Reaction Diagnostics for Ion species: " << pVacuumConfig->IonSpecies.Name << std::endl;
		}
	}

	out.flush();

	if ( pVacuumConfig->OutputFile != "" )
	{
		delete p_out;
	}


}

