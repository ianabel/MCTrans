
#include "NetCDFIO.hpp"
#include "MirrorPlasma.hpp"

// Code for NetCDF interface
//

using namespace netCDF;

NetCDFIO::NetCDFIO() 
{
}

void NetCDFIO::Open( const std::string &file )
{
	filename = file;
	data_file.open( file, netCDF::NcFile::FileMode::replace );
	TimeDim = data_file.addDim( "t" );
	TimeVar = data_file.addVar( "t", netCDF::NcDouble(), TimeDim );
	TimeVar.putAtt( "description", "Time since start of simulation" );
	TimeVar.putAtt( "units", "s" );
	TimeVar.putVar( {0}, 0.0 );
}

void NetCDFIO::Close()
{
	filename = "";
	data_file.close();
}

NetCDFIO::~NetCDFIO()
{
	if ( filename != "" )
		Close();
}

void NetCDFIO::AddScalarVariable( std::string name, std::string description, std::string units, double value )
{
	NcVar newvar = data_file.addVar( name, netCDF::NcDouble() );
	newvar.putAtt( "description", description );
	if ( units != "" )
		newvar.putAtt( "units", units );
	double tmp = value;
	newvar.putVar( &tmp );
}

void NetCDFIO::AddTextVariable( std::string name, std::string description, std::string units, std::string text )
{
	NcVar newvar = data_file.addVar( name, netCDF::NcString() );
	newvar.putAtt( "description", description );
	if ( units != "" )
		newvar.putAtt( "units", units );
	std::string value( text );
	const char *pStr = value.c_str();
	newvar.putVar( &pStr );
}

void NetCDFIO::AddTimeSeries( std::string name, std::string description, std::string units, double InitialValue )
{
	NcVar newvar = data_file.addVar( name, netCDF::NcDouble(), TimeDim );
	newvar.putAtt( "description", description );
	if ( units != "" )
		newvar.putAtt( "units", units );
	newvar.putVar( {0}, InitialValue );
}

size_t NetCDFIO::AddTimeSlice( double T )
{
	size_t next = TimeDim.getSize();
	std::vector<size_t> v = {next};
	TimeVar.putVar( v, T );
	return next;
}

void NetCDFIO::AppendToTimeSeries( std::string const& name, double value, size_t tIndex )
{
	NcVar variable = data_file.getVar( name );
	std::vector<size_t> v = {tIndex};
	variable.putVar( v, value );
}

// MirrorPlasma routines that use NetCDFIO
//

void MirrorPlasma::InitialiseNetCDF()
{
	if ( !isTimeDependent )
		return;
		
	if ( pVacuumConfig->NetcdfOutputFile == "" )
		return;

	nc_output.Open( pVacuumConfig->NetcdfOutputFile );
	nc_output.AddScalarVariable( "R_min","Innermost plasma radius", "m", pVacuumConfig->PlasmaInnerRadius() );
	nc_output.AddScalarVariable( "R_max","Outermost plasma radius", "m", pVacuumConfig->PlasmaOuterRadius() );
	nc_output.AddScalarVariable( "R_plasma","Plasma radius on centreline", "m", pVacuumConfig->PlasmaCentralRadius() );
	nc_output.AddScalarVariable( "IonDensity","Density of bulk ion species in the central cell", "10^20 m^-3", IonDensity );
	nc_output.AddScalarVariable( "ElectronDensity","Density of electrons in the central cell", "10^20 m^-3", ElectronDensity );
	nc_output.AddScalarVariable( "MirrorRatio","Ratio of Minimum to Maximum B along a field line", "", pVacuumConfig->MirrorRatio );

	// Time Dependent Variables
	nc_output.AddTimeSeries( "Voltage", "Voltage drop across the plasma","V", pVacuumConfig->ImposedVoltage );
	nc_output.AddTimeSeries( "MachNumber", "Plasma velocity divided by Sqrt(T_e/m_i)", "", MachNumber );
	nc_output.AddTimeSeries( "IonTemperature", "Temperature of the bulk ion species", "keV", IonTemperature );
	nc_output.AddTimeSeries( "ElectronTemperature", "Temperature of the bulk ion species", "keV", ElectronTemperature );
	nc_output.AddTimeSeries( "Current", "Radial current through the plasma","A", RadialCurrent() );


}

void MirrorPlasma::WriteTimeslice( double tNew )
{
	if ( !isTimeDependent )
		return;

	if ( pVacuumConfig->NetcdfOutputFile == "" )
		return;

	int tIndex = nc_output.AddTimeSlice( tNew );
	nc_output.AppendToTimeSeries( "Voltage", pVacuumConfig->ImposedVoltage, tIndex );
	nc_output.AppendToTimeSeries( "MachNumber", MachNumber, tIndex );
	nc_output.AppendToTimeSeries( "IonTemperature", IonTemperature, tIndex );
	nc_output.AppendToTimeSeries( "ElectronTemperature", ElectronTemperature, tIndex );
	nc_output.AppendToTimeSeries( "Current", RadialCurrent(), tIndex );
}

void MirrorPlasma::FinaliseNetCDF()
{
	if ( !isTimeDependent )
		return;

	if ( pVacuumConfig->NetcdfOutputFile == "" )
		return;

	nc_output.Close();
}


void MirrorPlasma::ReadVoltageFile( std::string const& fName )
{
	NcFile voltage_file( fName, NcFile::FileMode::read );
	if ( voltage_file.isNull() )
		throw std::runtime_error( "Could Not Open Voltage File " + fName + " as a NetCDF file" );

	NcVar time = voltage_file.getVar( "Time" );
	NcDim time_dim = voltage_file.getDim( "Time" );
	NcVar volt = voltage_file.getVar( "Voltage" );
	size_t nEntries = time_dim.getSize();
	std::vector<double> time_data( nEntries ),volt_data( nEntries );


	time.getVar( time_data.data() );
	volt.getVar( volt_data.data() );

#ifdef DEBUG
	std::cerr << "Last time in voltage file is " << time_data.back() << " with voltage " << volt_data.back() << std::endl;
#endif

	constexpr double epsilon = 1e-10; // extend the data just beyond the simulation end time to cope with floating point
	time_data.emplace_back( time_data.back() + epsilon );
	volt_data.emplace_back( volt_data.back() ); // Extend with final value

	VoltageFunction = std::make_unique<interpolant>( std::move( time_data ), std::move( volt_data ) );
}

void MirrorPlasma::WriteNetCDFReport(std::map<std::string, double>* parameterMap, int currentRun, int totalRuns)
{

	/*
	 *
	 */

	if ( isTimeDependent )
		return;

	if ( pVacuumConfig->NetcdfOutputFile == "" )
		return;

	std::string outputFileName = pVacuumConfig->NetcdfOutputFile;

	if ( totalRuns > 1 )
		outputFileName = pVacuumConfig->NetcdfOutputFile.insert( std::min(pVacuumConfig->NetcdfOutputFile.find_last_of("."), pVacuumConfig->NetcdfOutputFile.size()) , "_" + std::to_string(currentRun + 1));

	nc_output.Open( outputFileName );

	if( totalRuns > 1 && parameterMap != nullptr )
	{
		nc_output.AddScalarVariable( "BatchRunIndex","","",currentRun + 1.0 );
		for (auto const &pair: *parameterMap) 
		{
			nc_output.AddScalarVariable( pair.first, "Input Parameter","",pair.second );
		}
	}

	nc_output.AddTextVariable( "FuelName", "Fuel ion species", "", pVacuumConfig->IonSpecies.Name );

	nc_output.AddScalarVariable( "ZEff", "Effective ion charge", "", Zeff );
	nc_output.AddScalarVariable( "IonTemperature", "Temperature of the bulk ion species", "keV", IonTemperature );
	nc_output.AddScalarVariable( "ElectronTemperature", "Temperature of the bulk ion species", "keV", ElectronTemperature );
	nc_output.AddScalarVariable( "IonDensity","Density of bulk ion species in the central cell", "10^20 m^-3", IonDensity );
	nc_output.AddScalarVariable( "ElectronDensity","Density of electrons in the central cell", "10^20 m^-3", ElectronDensity );

	nc_output.AddScalarVariable( "R_min","Innermost plasma radius", "m", pVacuumConfig->PlasmaInnerRadius() );
	nc_output.AddScalarVariable( "R_max","Outermost plasma radius", "m", pVacuumConfig->PlasmaOuterRadius() );
	nc_output.AddScalarVariable( "R_plasma","Plasma mid-radius on midplane", "m", pVacuumConfig->PlasmaCentralRadius() );

	nc_output.AddScalarVariable( "R_wall", "Radius of central cell (first wall)", "m", pVacuumConfig->WallRadius );

	nc_output.AddScalarVariable( "R_throat", "Outer Radius of the plasma at the mirror throat", "m", (pVacuumConfig->PlasmaColumnWidth + pVacuumConfig->AxialGapDistance) /::sqrt( pVacuumConfig->MirrorRatio ) );

	nc_output.AddScalarVariable( "L_plasma", "Axial Length of the plasma", "m", pVacuumConfig->PlasmaLength );
	nc_output.AddScalarVariable( "V_plasma", "Plasma Volume", "m^3", pVacuumConfig->PlasmaVolume() );

	nc_output.AddScalarVariable( "B_midplane", "Magnetic Field Strength in the Central Cell", "T", pVacuumConfig->CentralCellFieldStrength );
	nc_output.AddScalarVariable( "B_throat", "Magnetic Field Strength in the mirror throat", "T", pVacuumConfig->CentralCellFieldStrength * pVacuumConfig->MirrorRatio );
	
	nc_output.AddScalarVariable( "RhoIon", "Ion Larmor Radius in the central cell", "m", IonLarmorRadius() );
	nc_output.AddScalarVariable( "a", "Typical plasma scale length", "m", pVacuumConfig->PlasmaColumnWidth/2.0 );
	nc_output.AddScalarVariable( "RhoStar", "Normalised larmor radius - rho_i/a", "", pVacuumConfig->PlasmaColumnWidth / ( 2.0 * IonLarmorRadius() ) );

	nc_output.AddScalarVariable( "AuxiliaryHeating", "Auxiliary electron heating.", "W", pVacuumConfig->AuxiliaryHeating * 1e6 );

	if ( pVacuumConfig->AlphaHeating ) {
		nc_output.AddScalarVariable( "AlphaHeating", "Electron heating from alpha particles", "MW", AlphaHeating()*pVacuumConfig->PlasmaVolume()*1e6 )  ;
		nc_output.AddScalarVariable( "AlphaPromptLossRate", "Prompt loss rate of energetic alpha particles", "s^-1", PromptAlphaLossFraction() * AlphaProductionRate() * pVacuumConfig->PlasmaVolume() );
		nc_output.AddScalarVariable( "AlphaPromptHeating", "Energy loss rate from prompt alpha particles", "W", AlphaPromptLosses() * pVacuumConfig->PlasmaVolume() * 1e6 );
		nc_output.AddScalarVariable( "AlphaSlowingDown", "Alpha particle slowing down time", "s", SlowingDownTime() );
	}

	nc_output.AddScalarVariable( "Mach", "Operating Mach Number", "", MachNumber );
	nc_output.AddScalarVariable( "AlfvenMach", "Alfven Mach Number", "", AlfvenMachNumber() );
	double v = SoundSpeed() * MachNumber;
	nc_output.AddScalarVariable( "Velocity","Peak rotational velocity", "m/s", v );
	double Rmid = pVacuumConfig->PlasmaCentralRadius();
	nc_output.AddScalarVariable( "AngularVelocity","Angular veloticy at the plasma centre", "s^-1", v/Rmid );

	nc_output.AddScalarVariable( "ViscousHeating", "Viscous Heating", "W", ViscousHeating() * pVacuumConfig->PlasmaVolume() );

	if ( Zeff > 0.0 )
	{
		nc_output.AddScalarVariable( "BremsstrahlungLosses", "Bremsstrahlung losses","W",BremsstrahlungLosses()*pVacuumConfig->PlasmaVolume() );
		nc_output.AddScalarVariable( "CyclotronLosses","Cyclotron emission losses","W", CyclotronLosses()*pVacuumConfig->PlasmaVolume() );
	}

	if ( pVacuumConfig->IncludeCXLosses ) {
		nc_output.AddScalarVariable( "CXLosses", "Heat loss due to charge exchange with neutrals", "W", CXHeatLosses()*pVacuumConfig->PlasmaVolume() );
	}
		
	nc_output.AddScalarVariable( "Voltage", "Total potential difference across the plasma","V",ElectricPotential() );

	double JRadial = ::fabs( RadialCurrent() );
	nc_output.AddScalarVariable( "JRadial", "Radial Current Drawn from Power Supply", "A", JRadial );
	nc_output.AddScalarVariable( "RotationPower", "Power Required (at the plasma) to support rotation","W", ElectricPotential() * JRadial );
	double omega = v/Rmid;
	nc_output.AddScalarVariable( "ViscousTorque","Power Loss from viscous torque", "W", ViscousTorque()*omega*pVacuumConfig->PlasmaVolume() );
	nc_output.AddScalarVariable( "ParallelRotationLosses","Power Loss from parallel losses", "W", ParallelAngularMomentumLossRate()*omega*pVacuumConfig->PlasmaVolume() );

	if (  pVacuumConfig->IncludeCXLosses ) {
		nc_output.AddScalarVariable( "CXRotationLosses","Rotational power loss from charge exchange", "W", CXMomentumLosses()*omega*pVacuumConfig->PlasmaVolume() );
	} 

	double KineticStoredEnergy = KineticEnergy();
	double ThermalStoredEnergy = ThermalEnergy();
	nc_output.AddScalarVariable( "KineticEnergy", "Plasma kinetic energy", "J", KineticStoredEnergy );
	nc_output.AddScalarVariable( "ThermalEnergy", "Thermal energy", "J", ThermalStoredEnergy );
	nc_output.AddScalarVariable( "StoredEnergy", "Total stored energy", "J", KineticStoredEnergy + ThermalStoredEnergy );

	nc_output.AddScalarVariable( "TauE", "Energy Confinement Time", "s", EnergyConfinementTime() );
	
	double ParallelConfinementTime = ( 1.5  * ElectronDensity * ElectronTemperature + 1.5 * IonDensity * IonTemperature ) * (  ReferenceDensity * ReferenceTemperature )/( ParallelElectronHeatLoss() + ParallelIonHeatLoss() );
	double PerpConfinementTime = ( 1.5  * IonDensity * IonTemperature * ReferenceDensity * ReferenceTemperature )/ClassicalIonHeatLoss();
	double CXConfinementTime = IonDensity*ReferenceDensity / CXLossRate();

	nc_output.AddScalarVariable( "ParallelConfinementTime", "Confinement time from parallel losses", "s", ParallelConfinementTime );
	nc_output.AddScalarVariable( "PerpConfinementTime", "Confinement time from perpendicular losses", "s", PerpConfinementTime );
	if ( pVacuumConfig->IncludeCXLosses ) {
		nc_output.AddScalarVariable( "CXConfinementTime", "Confinement time from charge-exchange losses", "s", CXConfinementTime );
	}

	double ParticleConfinementTime = ( IonDensity * ReferenceDensity ) /(  ParallelIonParticleLoss() + ClassicalIonParticleLosses() + CXLossRate() );

	nc_output.AddScalarVariable( "ParticleConfinementTime", "Particle (Ion) Confinement Time ~= ", "s", ParticleConfinementTime );

	nc_output.AddScalarVariable( "EquilibrationTime", "Ion-Electron Temperature Equilibration Time", "s", CollisionalTemperatureEquilibrationTime() );


	double ElectronLossRate = ParallelElectronParticleLoss() + ClassicalElectronParticleLosses();
	nc_output.AddScalarVariable( "ParticleLossRate","Total plasma loss rate ","electrons/s", ElectronLossRate*pVacuumConfig->PlasmaVolume() );

	nc_output.AddScalarVariable( "NeutralDensity","Density of neutral particles", "m^-3", NeutralDensity * ReferenceDensity );

	if ( !FixedNeutralDensity ) {
		ComputeSteadyStateNeutrals();
	}

	if ( pVacuumConfig->IncludeCXLosses ) {

		double CXR = CXLossRate() * pVacuumConfig->PlasmaVolume();
		nc_output.AddScalarVariable( "CXRate","Rate of charge-exchange reactions","s^-1", CXR );
		nc_output.AddScalarVariable( "CXHeatLoss","Thermal energy loss from charge-exchange","W", CXR * IonTemperature * ReferenceTemperature );

	} 

	double Resistance = ElectricPotential() / JRadial;
	double Capacitance = 2.0*KineticStoredEnergy / ( ElectricPotential() * ElectricPotential() );

	nc_output.AddScalarVariable( "Resistance","Electrical resistance of the plasma","Ohms", Resistance );
	nc_output.AddScalarVariable( "Capacitance","Electrical capacitance of the plasma","Farads", Capacitance );

	nc_output.AddScalarVariable( "Beta","Plasma beta","%",Beta()*100 );
	nc_output.AddScalarVariable( "IonCollisionality","Ratio of ion collision frequency to ion bounce time","",NuStar() );
	nc_output.AddScalarVariable( "HallParameter","Ratio of ion collision frequency to cyclotron frequency ( Omega_i * tau_ii)","",IonCyclotronFrequency()*IonCollisionTime() );
	
	double TripleProduct = IonDensity * ReferenceDensity * IonTemperature * EnergyConfinementTime();
	nc_output.AddScalarVariable( "TripleProduct", "n_i T_i τ_E", "keV s/m^3", TripleProduct );
	
	if ( pVacuumConfig->ReportThrust ) {
			nc_output.AddScalarVariable( "IonThrust","Thrust provided by fuel ions lost from the central cell","", ParallelIonThrust() );
			if ( pVacuumConfig->AlphaHeating )
				nc_output.AddScalarVariable( "AlphaThrust","Thrust provided by alpha particles promptly lost from the central cell","", PromptAlphaThrust() );
	}

	if ( pVacuumConfig->ReportNuclearDiagnostics ) {
		if ( pVacuumConfig->IonSpecies.Name == "Deuterium" ) {
			nc_output.AddScalarVariable( "DDNeutrons", "D/D Neutrons (2.45 MeV) per second","n/s", DDNeutronRate() ); 
			double Yield = ( 14.1/3.52 + 1.0 ) * FusionAlphaPowerDensity()*pVacuumConfig->PlasmaVolume();
			nc_output.AddScalarVariable( "DTEquivalentQ","Effective D/T Figure of Merit, assuming an identical plasma. Scientific D-T equivalent Q.","",Yield/( ElectricPotential()*JRadial ) );
		} else if ( pVacuumConfig->IonSpecies.Name == "Deuterium/Tritium Fuel" ) {

			double FusionAlphaPower = FusionAlphaPowerDensity()*pVacuumConfig->PlasmaVolume();
			double FusionNeutronPower = ( 14.1/3.52 ) * FusionAlphaPower;
			nc_output.AddScalarVariable( "FusionAlphaPower", "Fusion Power Output as Alphas","W",FusionAlphaPower*1e6 );
			nc_output.AddScalarVariable( "FusionNeutronPower", "Fusion Power Output as Alphas","W",FusionNeutronPower*1e6 );
			nc_output.AddScalarVariable( "NeutronWallLoading","Power impinging on the first wall as energetic neutrons","W/m^2",NeutronWallLoading()*1e6 );
			nc_output.AddScalarVariable( "ThermalOutputPower","Total Thermal Power Output","W", ThermalPowerOutput() *1e6 );

			double TotalOutputPower1 = FusionNeutronPower + FusionAlphaPower;
			double TotalOutputPower2 = ThermalPowerOutput();
			double TotalInputPower1 = ElectricPotential() * JRadial / 1e6; // Because the output is in MW
			double RotationDriveEta = 0.75;
			double TotalInputPower2 = ElectricPotential() * JRadial / ( RotationDriveEta * 1e6 ); // ditto

			nc_output.AddScalarVariable( "ScientificQ","Scientific Q-factor, not including tritium breeding, or engineering efficiencies","", TotalOutputPower1/TotalInputPower1 );
		} else {
			std::cerr << "Nuclear Reaction Diagnostics requested but no diagnostics exist for ion species: " << pVacuumConfig->IonSpecies.Name << std::endl;
		}
	}

	nc_output.Close();

}
