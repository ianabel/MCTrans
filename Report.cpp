
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

void MirrorPlasma::PrintReport(std::map<std::string, double>* parameterMap, int currentRun, int totalRuns)
{
	std::ostream *p_out;
	if ( OutputFile == "" )
		p_out = &std::cout;
	else if(totalRuns == 1)
	{
		p_out = new std::fstream( OutputFile, std::ios_base::out | std::ios_base::trunc );
	}
	else
	{
		std::string outputFileName = OutputFile.insert( std::min(OutputFile.find_last_of("."), OutputFile.size()) , "_" + std::to_string(currentRun + 1));
		p_out = new std::fstream( outputFileName, std::ios_base::out | std::ios_base::trunc );
	}

	std::ostream &out = *p_out;

	out.precision( 3 );

	if( totalRuns > 1 && parameterMap)
	{
		out << "Batch run " << currentRun+1 << " of " << totalRuns << std::endl;
		out << std::endl;
		out << "Input parameters: " << std::endl;
		for (auto const &pair: *parameterMap) 
			out << pair.first << ": " << pair.second << std::endl;
		out << std::endl;
	}

	out << "The plasma is made up of electrons, " << IonSpecies.Name; 
	if ( Zeff > 0.0 )
		out << " and trace, radiating, impurities" << std::endl;
	else 
		out << std::endl;

	out << "Electron Density is " << ElectronDensity*ReferenceDensity << " m^-3" << std::endl;
	out << "Ion Density is " << IonDensity * ReferenceDensity << " m^-3" << std::endl;
	out << std::endl;
	out << "Electron Temperature is "; PrintWithUnit( out, ElectronTemperature * 1e3,"eV"); out << std::endl;
	out << "Ion Temperature is "; PrintWithUnit( out, IonTemperature * 1e3,"eV"); out << std::endl;

	out << "Radius of central cell (first wall) is " << WallRadius << " m" << std::endl;

	out << "\t Inner radius of the plasma is " << AxialGapDistance << " m" << std::endl;
	out << "\t Outer radius of the plasma is " << PlasmaColumnWidth + AxialGapDistance << " m" << std::endl;

	out << "Outer Radius of the plasma at the mirror throat is " << (PlasmaColumnWidth + AxialGapDistance) /::sqrt( MirrorRatio ) << " m" << std::endl;
	out << std::endl;
	out << "Axial Length of the plasma is " << PlasmaLength << " m" << std::endl;
	out << "Plasma Volume is " << PlasmaVolume() << " m^3" << std::endl;
	out << "" << std::endl;
	out << "Magnetic Field in the Central Cell is " << CentralCellFieldStrength << " T" << std::endl;
	out << "Magnetic Field at the Mirror Throat is " << MirrorRatio * CentralCellFieldStrength << " T" << std::endl;
	out << "" << std::endl;
	out << "Ion Larmor Radius in the central cell is " << IonLarmorRadius() << " m " << std::endl;
	out << "Typical plasma scale lengths are ~ " << PlasmaColumnWidth/2.0 << " m  = " << PlasmaColumnWidth / ( 2.0 * IonLarmorRadius() ) << " rho_i" << std::endl;

	out << std::endl;

	if ( AuxiliaryHeating > 0 ) {
		out << "Auxiliary Heating of " << AuxiliaryHeating * 1e3  << " kW was included" << std::endl;
		out << "" << std::endl;
	} else {
		out << "No auxiliary heating was included in this calculation." << std::endl;
		out << std::endl;
	}

	if ( IncludeAlphaHeating ) {
		out << "Self-consistent Alpha Heating was included in this calculation. " << std::endl;

		out << "  Alphas provide " << AlphaHeating()*PlasmaVolume() << " MW of heating" << std::endl;
		out << "  Alpha particles are lost directly at a rate of " << PromptAlphaLossFraction() * AlphaProductionRate() * PlasmaVolume() << " /s" << std::endl;
		out << "  Prompt Alpha Losses towards the end plates give " << AlphaPromptLosses() * PlasmaVolume() << " MW of energy losses" << std::endl;

		double TauSD = SlowingDownTime();
		out << "  The alpha particle slowing-down time is " << TauSD << " s" << std::endl;
	}

	out << std::endl;
	out << "Operating Mach number is " << MachNumber << std::endl;
	out << "Alfven Mach number is " << AlfvenMachNumber() << std::endl;
	out << std::endl;
	double v = SoundSpeed() * MachNumber;
	out << "Velocity is " << v << " m/s" << std::endl;
	double Rmid = PlasmaCentralRadius();
	out << "Angular Velocity at the plasma centre (R = " << Rmid << " m) is " << v/Rmid << " /s" << std::endl;
	out << std::endl;

	out << "Viscous Heating is ";
	PrintWithUnit( out, ViscousHeating() * PlasmaVolume(),"W" ); out << std::endl;
	if ( Zeff > 0.0 )
	{
		out << "Bremsstrahlung losses are ";
		PrintWithUnit( out, BremsstrahlungLosses()*PlasmaVolume(), "W" ); out << std::endl;
		out << "Cyclotron emission losses are ";
		PrintWithUnit( out, CyclotronLosses()*PlasmaVolume(), "W" ); out << std::endl;
	}
	else
		out << "No radiation losses were included." << std::endl;

	if ( IncludeCXLosses ) {
		out << "Heat loss due to charge exchange with neutrals is ";
		PrintWithUnit( out, CXHeatLosses()*PlasmaVolume(), "W" ); out << std::endl;
	} else {
		out << "Losses due to charge exchange were not included" << std::endl;
	}

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
	out << "\t Power Loss from viscous torque  "; PrintWithUnit( out, ViscousTorque()*omega*PlasmaVolume(), "W" ); out << std::endl;
	out << "\t Power Loss from parallel loss   "; PrintWithUnit( out, ParallelAngularMomentumLossRate()*omega*PlasmaVolume(), "W" ); out << std::endl;
	if (  IncludeCXLosses ) {
		out << "\t Power Loss from charge exchange "; PrintWithUnit( out, CXMomentumLosses()*omega*PlasmaVolume(), "W" ); out << std::endl;
	} else {
		out << "\t Power Loss from charge exchange was not included. " << std::endl;
	}
	out << std::endl;

	double KineticStoredEnergy = KineticEnergy();
	double ThermalStoredEnergy = ThermalEnergy();
	out << "Total Stored Energy is "; PrintWithUnit( out, KineticStoredEnergy + ThermalStoredEnergy, "J" ); out << std::endl;
	out << " of which "; PrintWithUnit( out, KineticStoredEnergy, "J" ); out << " is kinetic energy in the rotation" << std::endl;
	out << "      and "; PrintWithUnit( out, ThermalStoredEnergy, "J" ); out << " is thermal energy of the plasma" << std::endl;


	out << "Energy Confinement Time is "; PrintWithUnit( out, EnergyConfinementTime(),"s" ); out << std::endl;
	out << "Of which" << std::endl;
	
	double ParallelConfinementTime = ( 1.5  * ElectronDensity * ElectronTemperature + 1.5 * IonDensity * IonTemperature ) * (  ReferenceDensity * ReferenceTemperature )/( ParallelElectronHeatLoss() + ParallelIonHeatLoss() );
	double PerpConfinementTime = ( 1.5  * IonDensity * IonTemperature * ReferenceDensity * ReferenceTemperature )/( ClassicalIonHeatLoss() + ClassicalElectronHeatLoss() );
	double CXConfinementTime = IonDensity*ReferenceDensity / CXLossRate();
	out << "\tConfinement time from parallel losses is\t"; PrintWithUnit( out, ParallelConfinementTime, "s" ); out << std::endl;
	out << "\tConfinement time from perpendicular losses is\t"; PrintWithUnit( out, PerpConfinementTime, "s" ); out << std::endl;
	if ( IncludeCXLosses ) {
		out << "\tConfinement time from charge-exchange losses is\t"; PrintWithUnit( out, CXConfinementTime, "s" ); out << std::endl;
	}
	out << std::endl;

	double ParticleConfinementTime = ( IonDensity * ReferenceDensity ) /(  ParallelIonParticleLoss() + ClassicalIonParticleLosses() + CXLossRate() );
	out << "Particle (Ion) Confinement Time ~= "; PrintWithUnit( out, ParticleConfinementTime, "s" ); out << std::endl;
	out << std::endl;

	out << "Ion-Electron Temperature Equilibration Time is "; PrintWithUnit( out, CollisionalTemperatureEquilibrationTime(),"s" ); out << std::endl;

	out << std::endl;
	double ElectronLossRate = ParallelElectronParticleLoss() + ClassicalElectronParticleLosses();
	out << "Plasma must be provided at a rate of " << ElectronLossRate*PlasmaVolume() << " electrons /s to maintain steady-state" << std::endl;

	out << std::endl;
	if ( FixedNeutralDensity ) {
		out << "A fixed neutral density of " << NeutralDensity * ReferenceDensity << " /m^3 has been assumed." << std::endl;
	} else {
		ComputeSteadyStateNeutrals();
		out << "Neutral gas must be provided at a rate of " << NeutralSource * PlasmaVolume() << " particles/s to refuel the plasma" << std::endl;
		out << "This leads to a steady-state neutral density of " << NeutralDensity * ReferenceDensity << "/m^3" << std::endl;
	}
	if ( IncludeCXLosses ) {
		double CXR = CXLossRate() * PlasmaVolume();
		out << "The level of charge-exchange losses due to the neutrals is " << CXR << " particles/s lost" << std::endl;
		out << "\t This is equivalent to a heat loss rate of "; PrintWithUnit( out, CXR * IonTemperature * ReferenceTemperature, "W" );out << std::endl;
	} else {
		out << "Steady-state charge exchange losses were not included." << std::endl;
	}


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
	// if ( NuStar() > 1.0 )
	//	out << "\t\t 
	out << "\t ρ* = " << ( 2.0 * IonLarmorRadius() ) / PlasmaColumnWidth << std::endl; 
	out << "\t Omega_i tau_ii = " << IonCyclotronFrequency()*IonCollisionTime() << std::endl; 

	

	out << std::endl;
	out << "Fusion Triple Product " << std::endl;
	double TripleProduct = IonDensity * ReferenceDensity * IonTemperature * EnergyConfinementTime();
	out << "\t n T τ = " << TripleProduct << " keV s /m^3" << std::endl;
	
	out << std::endl;

	if ( ReportThrust ) {
			out << "Ions lost from the central cell provide " << ParallelIonThrust() << " N of thrust" << std::endl;
			if ( IncludeAlphaHeating )
				out << "Prompt Alpha losses contribute " << PromptAlphaThrust() << " N of thrust" << std::endl;
	}


	if ( ReportNuclearDiagnostics ) {
		if ( IonSpecies.Name == "Deuterium" ) {
			out << " === Nuclear Reactions Assuming Deuterium Fuel === " << std::endl;
			out << " D/D Neutrons (2.45 MeV) per second: " << DDNeutronRate() << std::endl;
			// out << " D/T Neutrons (14.1 MeV) from Fusion-Produced-Tritium: " << DDEnergeticNeutronRate() << std::endl;
			out << std::endl;

			out << " === Effective D/T Figures of Merit, assuming an identical plasma === " << std::endl;
			double Yield = ( 14.1/3.52 + 1.0 ) * FusionAlphaPowerDensity()*PlasmaVolume();
			out << " Q_(Scientific DT equivalent) = " << Yield / ElectricPotential()*JRadial << std::endl;
		} else if ( IonSpecies.Name == "Deuterium/Tritium Fuel" ) {
			out << " === Reactor Output Assuming D/T Reactor Fuel === " << std::endl;

			double FusionAlphaPower = FusionAlphaPowerDensity()*PlasmaVolume();
			double FusionNeutronPower = ( 14.1/3.52 ) * FusionAlphaPower;
			out << "Fusion Power Output as Alphas:   " << FusionAlphaPower   << " MW" << std::endl;
			out << "                    as Neutrons: " << FusionNeutronPower << " MW" << std::endl;
			out << "Neutron Wall Loading is " << NeutronWallLoading() << " MW/m^2" << std::endl;
			out << "Total Thermal Power Output is " << ThermalPowerOutput() << " MW" << std::endl;

			out << " === Figures of Merit === " << std::endl;
			double TotalOutputPower1 = FusionNeutronPower + FusionAlphaPower;
			double TotalOutputPower2 = ThermalPowerOutput();
			double TotalInputPower1 = ElectricPotential() * JRadial / 1e6; // Because the output is in MW
			double RotationDriveEta = 0.75;
			double TotalInputPower2 = ElectricPotential() * JRadial / ( RotationDriveEta * 1e6 ); // ditto

			out << "\t Q_scientific (not including Tritium Breeding or efficiencies ) = " << TotalOutputPower1/TotalInputPower1 << std::endl;
			out << "\t Q_engineering (including Tritium Breeding at TBR of 1.0 and assuming efficiency " << RotationDriveEta << " for electric drive of rotation ) = " << TotalOutputPower2/TotalInputPower2 << std::endl;

		} else {
			out << "No Nuclear Reaction Diagnostics for Ion species: " << IonSpecies.Name << std::endl;
		}
	}

	out.flush();

	if ( OutputFile != "" )
	{
		delete p_out;
	}


}

