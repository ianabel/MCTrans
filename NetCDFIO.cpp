
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

NetCDFIO::~NetCDFIO()
{
}

void NetCDFIO::AddScalarVariable( std::string&& name, std::string&& description, std::string&& units, double value )
{
	NcVar newvar = data_file.addVar( name, netCDF::NcDouble() );
	newvar.putAtt( "description", description );
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
	nc_output.AddScalarVariable( "IonDensity","Density of bulk ion species in the central cell", "m^-3", IonDensity );
	nc_output.AddScalarVariable( "ElectronDensity","Density of electrons in the central cell", "m^-3", ElectronDensity );
	
}
