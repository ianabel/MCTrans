
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

}

void MirrorPlasma::WriteTimeslice()
{
	// time dependent (or will be soon)
	nc_output.AddScalarVariable( "Voltage", "Voltage drop across the plasma","V", pVacuumConfig->ImposedVoltage );
	nc_output.AddScalarVariable( "MachNumber", "Plasma velocity divided by Sqrt(T_e/m_i)", "", MachNumber );
	nc_output.AddScalarVariable( "IonTemperature", "Temperature of the bulk ion species", "keV", IonTemperature );
	nc_output.AddScalarVariable( "ElectronTemperature", "Temperature of the bulk ion species", "keV", ElectronTemperature );
	nc_output.AddScalarVariable( "Current", "Radial current through the plasma","A", RadialCurrent() );

}

void MirrorPlasma::FinaliseNetCDF()
{
	if ( pVacuumConfig->NetcdfOutputFile == "" )
		return;

	nc_output.Close();
}
