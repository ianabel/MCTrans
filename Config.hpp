
#include "Centrifugal.hpp"
#include "toml11/toml.hpp"

void LoadConfig( std::string const& configFile, Plasma &plasma, Configuration& conf )
{
	const toml::value config = toml::parse( configFile );

	if ( config.count( "plasma" ) != 1 )
		throw std::invalid_argument( "There must be one, and only one, [plasma] block in the configuration file." );

	const auto plasmaConfig = toml::find( config, "plasma" );

	if ( plasmaConfig.count( "IonSpecies" ) != 1 )
		throw std::invalid_argument( "Fuel must be specified once in the [plasma] block" );

	std::string FuelName = plasmaConfig.at( "IonSpecies" ).as_string();


	if ( FuelName == "Hydrogen" ) {
		plasma.Mu = 1.0;
		plasma.ZIon = 1.0;
		plasma.Fuel = "Hydrogen";
		conf.ReportNuclear = false;
		conf.IncludeAlphaHeating = false;
	} else if ( FuelName == "Deuterium" ) {
		plasma.Mu = 2.0;
		plasma.ZIon = 1.0;
		plasma.Fuel = "Deuterium";
		conf.ReportNuclear = false;
		conf.IncludeAlphaHeating = false;
	} else if ( FuelName == "DT Fuel" ) {
		plasma.Mu = 2.5;
		plasma.ZIon = 1.0;
		plasma.Fuel = "Deuterium/Tritium Fuel";
		conf.ReportNuclear = true;
		conf.IncludeAlphaHeating = true;
	}

	plasma.TiTe = toml::find<double>( plasmaConfig, "IonToElectronTemperatureRatio" );
	plasma.Zeff = toml::find<double>( plasmaConfig, "Zeff" );

	plasma.LumpedImpurity = false;

	plasma.ElectronDensity = toml::find<double>( plasmaConfig, "ElectronDensity" );

	plasma.IonDensity = plasma.ElectronDensity / plasma.ZIon;
	if ( FuelName == "DT Fuel" )
		plasma.IonDensity /= 2.0;
	

	conf.MachSolve = false;
	conf.TempSolve = false;
	if ( plasmaConfig.count( "ElectronTemperature" ) == 1 )
	{
		plasma.ElectronTemperature = toml::find<double>( plasmaConfig, "ElectronTemperature" );
		conf.MachSolve = true;
	}

	const auto mirrorConfig = toml::find<toml::table>( config, "configuration" );


	conf.CentralCellFieldStrength =  mirrorConfig.at( "CentralCellField" ).as_floating();

	if ( mirrorConfig.count( "MirrorRatio" ) == 1 ) {
		// set Mirror Ratio directly
		conf.MirrorRatio = mirrorConfig.at( "MirrorRatio" ).as_floating();

		if ( mirrorConfig.count( "ThroatField" ) == 1 ) {
			throw std::invalid_argument( toml::format_error( "[error] Cannot specify botth Miror Ratio and Throat Field",mirrorConfig.at( "MirrorRatio" )," mirror ratio defined here", mirrorConfig.at( "ThroatField" ), " Throat field here" ) );
		}
	} else if ( mirrorConfig.count( "ThroatField" ) == 1 ) {
		double MagFieldThroat = mirrorConfig.at( "ThroatField" ).as_floating();
		conf.MirrorRatio = MagFieldThroat / conf.CentralCellFieldStrength;
		if ( mirrorConfig.count( "MirrorRatio" ) == 1 ) {
			throw std::invalid_argument( toml::format_error( "[error] Cannot specify botth Miror Ratio and Throat Field",mirrorConfig.at( "MirrorRatio" )," mirror ratio defined here", mirrorConfig.at( "ThroatField" ), " Throat field here" ) );
		}
	} else {
		throw std::invalid_argument( "[error] Must specify either MirrorRatio or ThroatField" );
	}

	if ( mirrorConfig.count( "PlasmaRadiusMin" ) == 1  || mirrorConfig.count( "PlasmaRadiusMax" ) == 1 ) {
		if ( mirrorConfig.count( "PlasmaRadiusMin" ) == 0 )
		{
			throw std::invalid_argument( toml::format_error( "[error] When PlasmaRadiusMax is specified you must also set PlasmaRadiusMin",mirrorConfig.at( "PlasmaRadiusMax" )," max radius set here" ) );
		}

		if ( mirrorConfig.count( "PlasmaRadiusMax" ) == 0 )
		{
			throw std::invalid_argument( toml::format_error( "[error] When PlasmaRadiusMin is specified you must also set PlasmaRadiusMax",mirrorConfig.at( "PlasmaRadiusMin" )," min radius set here" ) );
		}
		// set Plasma Radii directly
		double PlasmaMinRadius = mirrorConfig.at( "PlasmaMinRadius" ).as_floating();
		double PlasmaMaxRadius = mirrorConfig.at( "PlasmaMaxRadius" ).as_floating();

		conf.AxialGapDistance = PlasmaMinRadius;
		conf.PlasmaColumnWidth = PlasmaMaxRadius - PlasmaMinRadius;


		if ( mirrorConfig.count( "AxialGapDistance" ) == 1 ) {
			throw std::invalid_argument( toml::format_error( "[error] When PlasmaRadiusMin / Max are specified you cannot set AxialGapDistance",mirrorConfig.at( "PlasmaRadiusMin" )," minimum radius set here", mirrorConfig.at( "AxialGapDistance" ), " AxialGapDistance found here" ) );
		} else if ( mirrorConfig.count( "PlasmaColumnWidth") == 1 ) {
			throw std::invalid_argument( toml::format_error( "[error] When PlasmaRadiusMin / Max are specified you cannot set AxialGapDistance",mirrorConfig.at( "PlasmaRadiusMin" )," minimum radius set here", mirrorConfig.at( "AxialGapDistance" ), " AxialGapDistance found here" ) );
		}
	} else if ( mirrorConfig.count( "AxialGapDistance" ) == 1  || mirrorConfig.count( "PlasmaColumnWidth" ) == 1 ) {
		if ( mirrorConfig.count( "AxialGapDistance" ) == 0 )
		{
			throw std::invalid_argument( toml::format_error( "[error] When PlasmaColumnWidth is specified you must also set AxialGapDistance",mirrorConfig.at( "PlasmaColumnWidth" )," width set here" ) );
		}

		if ( mirrorConfig.count( "PlasmaColumnWidth" ) == 0 )
		{
			throw std::invalid_argument( toml::format_error( "[error] When AxialGapDistance is specified you must also set PlasmaColumnWidth",mirrorConfig.at( "AxialGapDistance" )," axial gap set here" ) );
		}
		// set parameters directly

		conf.AxialGapDistance = mirrorConfig.at( "AxialGapDistance" ).as_floating();
		conf.PlasmaColumnWidth = mirrorConfig.at( "PlasmaColumnWidth" ).as_floating();



		if ( mirrorConfig.count( "PlasmaMinRadius" ) == 1 ) {
			throw std::invalid_argument( toml::format_error( "[error] When AxialGapDistance and PlasmaColumnWidth are specified you cannot set PlasmaMinRadius",mirrorConfig.at( "PlasmaRadiusMin" )," minimum radius set here", mirrorConfig.at( "AxialGapDistance" ), " AxialGapDistance found here" ) );
		} else if ( mirrorConfig.count( "PlasmaMaxRadius") == 1 ) {
			throw std::invalid_argument( toml::format_error( "[error] When AxialGapDistance and PlasmaColumnWidth are specified you cannot set PlasmaMaxRadius",mirrorConfig.at( "PlasmaRadiusMax" )," maximum radius set here", mirrorConfig.at( "AxialGapDistance" ), " AxialGapDistance found here" ) );
		}
	} else {
		throw std::invalid_argument( "[error] Must specify plasma annulus with either PlasmaRadiusMin / PlasmaRadiusMax or AxialGapDistance and PlasmaColumnWidth" );
	}

	if ( mirrorConfig.count( "Voltage" ) != 0 ) {
		if ( conf.MachSolve ) {
			throw std::invalid_argument( "[error] Cannot specify both temperature and voltage." );
		} else {
			conf.Voltage = mirrorConfig.at( "Voltage" ).as_floating();
			conf.TempSolve = true;
		}
	}

	if ( mirrorConfig.count( "Voltage" ) == 0 && plasmaConfig.count( "ElectronTemperature" ) == 0 ) {
		throw std::invalid_argument( "[error] Must specify at least one of ElectronTemperature or Voltage." );
	}

	conf.WallRadius = mirrorConfig.at( "WallRadius" ).as_floating();

	conf.PlasmaLength = mirrorConfig.at( "PlasmaLength" ).as_floating();

	if ( mirrorConfig.count( "AuxiliaryHeating" ) == 1 )
		conf.AuxiliaryHeating = mirrorConfig.at( "AuxiliaryHeating" ).as_floating();
	else
		conf.AuxiliaryHeating = 0.0;

	// In case overridden
	if ( mirrorConfig.count( "ReportNuclearDiagnostics" ) == 1 )
		conf.ReportNuclear = mirrorConfig.at( "ReportNuclearDiagnostics" ).as_boolean();

	if ( mirrorConfig.count( "IncludeAlphaHeating" ) == 1 )
		conf.IncludeAlphaHeating = mirrorConfig.at( "IncludeAlphaHeating" ).as_boolean();

	if ( mirrorConfig.count( "NeutralDensity" ) == 1 ) {
		plasma.NeutralSource = 0;
		plasma.NeutralDensity = mirrorConfig.at( "NeutralDensity" ).as_floating();
	} else {
		plasma.NeutralSource = 0;
		plasma.NeutralDensity = 0;
	}

}
