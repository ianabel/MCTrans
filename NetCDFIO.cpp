
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
}

void NetCDFIO::AddScalarVariable( std::string&& name, std::string&& description, std::string&& units, double value )
{
	NcVar newvar = data_file.addVar( name, netCDF::NcDouble() );
	newvar.putAtt( "description", description );
	if ( units != "" )
		newvar.putAtt( "units", units );
	double tmp = value;
	newvar.putVar( &tmp );
}

void NetCDFIO::AddTimeSeries( std::string&& name, std::string&& description, std::string&& units, double InitialValue )
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
