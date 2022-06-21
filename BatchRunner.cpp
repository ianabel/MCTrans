#include <toml.hpp>
#include <iostream>

#include "BatchRunner.hpp"
#include "MirrorPlasma.hpp"
#include "Config.hpp"

// Type names used in the recursive cartesian product
using vecMap = std::vector<std::map<std::string, double>>;
using vecPairVS = std::vector<std::pair<std::vector<double>*,std::string>>;
using mapSD = std::map<std::string, double>;

BatchRunner::BatchRunner(std::string const& batchFile)
{
	const auto batchConfig = toml::parse( batchFile );
	// Algorithm Parameters
	if ( batchConfig.count( "algorithm" ) != 0 )
	{
		const auto algConfig = toml::find<toml::table>( batchConfig, "algorithm" );

		readParameterFromFile(algConfig, "ParallelFudgeFactor", ParallelFudgeFactorVals, false, 1.0);
		readParameterFromFile(algConfig, "PerpFudgeFactor", PerpFudgeFactorVals, false, 1.0);
		readParameterFromFile(algConfig, "InitialTemp", InitialTempVals, false, 0.1);
		readParameterFromFile(algConfig, "InitialMach", InitialMachVals, false, 4.0);
		readParameterFromFile(algConfig, "RateThreshold", RateThresholdVals, false, 1e-4);
		readParameterFromFile(algConfig, "SundialsAbsTol", SundialAbsTolVals, false, 1e-9);
		readParameterFromFile(algConfig, "SundialsRelTol", SundialsRelTolVals, false, 1e-7);
		
		if ( algConfig.count( "UseAmbipolarPhi" ) == 1 )
			AmbipolarPhi = algConfig.at( "UseAmbipolarPhi" ).as_boolean();
		else
			AmbipolarPhi = true;

		if ( algConfig.count( "AsciiOutputFile" ) == 1 )
		{
			OutputFile = algConfig.at( "AsciiOutputFile" ).as_string();
		} else {
			OutputFile = "";
		}

		if ( algConfig.count( "NetcdfOutput" ) == 1 )
		{
			NetcdfOutputFile = algConfig.at( "NetcdfOutput" ).as_string();
		} else {
			NetcdfOutputFile = "";
		}

		if ( algConfig.count( "IncludeChargeExchangeLosses" ) == 1 )
			IncludeCXLosses = algConfig.at( "IncludeChargeExchangeLosses" ).as_boolean();
		else
			IncludeCXLosses = false;

		if ( algConfig.count( "UseCollisionalFluxes" ) == 1 )
			Collisional = algConfig.at( "UseCollisionalFluxes" ).as_boolean();
		else
			Collisional = false;

#ifdef DEBUG
		if ( ( algConfig.count( "InitialTemp" ) == 1 ) )
			std::cerr << "Initial Temperature for Temperature Solve set from config file to " << InitialTempVals.size() << " value(s) between " << InitialTempVals[0] << " and " << InitialTempVals.back() << std::endl;
		else
			std::cerr << "Initial Temperature for Temperature Solve set to the default of " << InitialTempVals[0] << std::endl;

		if ( algConfig.count( "InitialMach" ) == 1 )
			std::cerr << "Initial Mach Number for fixed-temperature solve set from config file to " << InitialMachVals.size() << " value(s) between " << InitialMachVals[0] << " and " << InitialMachVals.back() << std::endl;
		else 
			std::cerr << "Initial Mach Number for fixed-temperature solve set to the default of " << InitialMachVals[0] << std::endl;
#endif
	}
	else
	{
#ifdef DEBUG
		std::cerr << "No [algorithm] section, using default values for internal knobs." << std::endl;
#endif
		ParallelFudgeFactorVals.push_back(1.0);
		PerpFudgeFactorVals.push_back(1.0);
		InitialTempVals.push_back(0.1);
		InitialMachVals.push_back(4.0);
		SundialAbsTolVals.push_back(1e-7);
		SundialsRelTolVals.push_back(1e-7);
		RateThresholdVals.push_back(1e-4);
		ptrsAndNamesToVectors.push_back(std::make_pair(&ParallelFudgeFactorVals, "ParallelFudgeFactor"));
		ptrsAndNamesToVectors.push_back(std::make_pair(&PerpFudgeFactorVals, "PerpFudgeFactor"));
		ptrsAndNamesToVectors.push_back(std::make_pair(&InitialTempVals, "InitialTemp"));
		ptrsAndNamesToVectors.push_back(std::make_pair(&InitialMachVals, "InitialMach"));
		ptrsAndNamesToVectors.push_back(std::make_pair(&SundialAbsTolVals, "SundialsAbsTol"));
		ptrsAndNamesToVectors.push_back(std::make_pair(&SundialAbsTolVals, "SundialsRelTol"));
		ptrsAndNamesToVectors.push_back(std::make_pair(&RateThresholdVals, "RateThreshold"));

		AmbipolarPhi = true;
		Collisional = false;
		IncludeCXLosses = false;
		OutputFile  = "";
		NetcdfOutputFile = "";
	}

	const auto batch = toml::find<toml::value>( batchConfig, "configuration" );
	// Fuel name
	if ( batch.count( "IonSpecies" ) != 1 )
		throw std::invalid_argument( "[error] Fuel must be specified once in the [configuration] block" );
	FuelName = batch.at( "IonSpecies" ).as_string();

	// Report Thrust
	if ( batch.count( "ReportThrust" ) == 1 )
		ReportThrust = batch.at( "ReportThrust" ).as_boolean();
	else
		ReportThrust = false;

	// Central Field Strength
	readParameterFromFile(batch, "CentralCellField", CentralCellFieldStrengthVals);
	
	// Mirror Ratio
	if ( batch.count( "MirrorRatio" ) == 1 ) 
	{
		readParameterFromFile(batch, "MirrorRatio", MirrorRatioVals, false);
		if ( batch.count( "ThroatField" ) == 1 ) {
			throw std::invalid_argument( toml::format_error( "[error] Cannot specify botth Miror Ratio and Throat Field",batch.at( "MirrorRatio" )," mirror ratio defined here", batch.at( "ThroatField" ), " Throat field here" ) );
		}
	}
	// Magnetic Throat Field
	else if ( batch.count( "ThroatField" ) == 1 ) 
	{
		readParameterFromFile(batch, "ThroatField", MagFieldThroatVals, false);
		if ( batch.count( "MirrorRatio" ) == 1 ) {
			throw std::invalid_argument( toml::format_error( "[error] Cannot specify botth Miror Ratio and Throat Field",batch.at( "MirrorRatio" )," mirror ratio defined here", batch.at( "ThroatField" ), " Throat field here" ) );
		}
	} 
	else {
		throw std::invalid_argument( "[error] Must specify either MirrorRatio or ThroatField" );
	}

	// Plasma Radius Max and Min
	if ( batch.count( "PlasmaRadiusMin" ) == 1  || batch.count( "PlasmaRadiusMax" ) == 1 ) {
		if ( batch.count( "PlasmaRadiusMin" ) == 0 )
			throw std::invalid_argument( toml::format_error( "[error] When PlasmaRadiusMax is specified you must also set PlasmaRadiusMin",batch.at( "PlasmaRadiusMax" )," max radius set here" ) );
		if ( batch.count( "PlasmaRadiusMax" ) == 0 )
			throw std::invalid_argument( toml::format_error( "[error] When PlasmaRadiusMin is specified you must also set PlasmaRadiusMax",batch.at( "PlasmaRadiusMin" )," min radius set here" ) );
		
		readParameterFromFile(batch, "PlasmaRadiusMin", AxialGapDistanceVals, false);
		readParameterFromFile(batch, "PlasmaRadiusMax", PlasmaMaxRadiusVals, false);

		if ( batch.count( "AxialGapDistance" ) == 1 ) {
			throw std::invalid_argument( toml::format_error( "[error] When PlasmaRadiusMin / Max are specified you cannot set AxialGapDistance",batch.at( "PlasmaRadiusMin" )," minimum radius set here", batch.at( "AxialGapDistance" ), " AxialGapDistance found here" ) );
		} else if ( batch.count( "PlasmaColumnWidth" ) == 1 ) {
			throw std::invalid_argument( toml::format_error( "[error] When PlasmaRadiusMin / Max are specified you cannot set AxialGapDistance",batch.at( "PlasmaRadiusMin" )," minimum radius set here", batch.at( "AxialGapDistance" ), " AxialGapDistance found here" ) );
		}
	}

	// Axial Gap Distance and Plasma Column Width
	else if ( batch.count( "AxialGapDistance" ) == 1  || batch.count( "PlasmaColumnWidth" ) == 1 ) {
		if ( batch.count( "AxialGapDistance" ) == 0 )
			throw std::invalid_argument( toml::format_error( "[error] When PlasmaColumnWidth is specified you must also set AxialGapDistance",batch.at( "PlasmaColumnWidth" )," width set here" ) );
		if ( batch.count( "PlasmaColumnWidth" ) == 0 )
			throw std::invalid_argument( toml::format_error( "[error] When AxialGapDistance is specified you must also set PlasmaColumnWidth",batch.at( "AxialGapDistance" )," axial gap set here" ) );
		
		readParameterFromFile(batch, "AxialGapDistance", AxialGapDistanceVals, false);
		readParameterFromFile(batch,"PlasmaColumnWidth", PlasmaColumnWidthVals, false);

		if ( batch.count( "PlasmaMinRadius" ) == 1 ) {
			throw std::invalid_argument( toml::format_error( "[error] When AxialGapDistance and PlasmaColumnWidth are specified you cannot set PlasmaMinRadius",batch.at( "PlasmaRadiusMin" )," minimum radius set here", batch.at( "AxialGapDistance" ), " AxialGapDistance found here" ) );
		} else if ( batch.count( "PlasmaMaxRadius" ) == 1 ) {
			throw std::invalid_argument( toml::format_error( "[error] When AxialGapDistance and PlasmaColumnWidth are specified you cannot set PlasmaMaxRadius",batch.at( "PlasmaRadiusMax" )," maximum radius set here", batch.at( "AxialGapDistance" ), " AxialGapDistance found here" ) );
		}
	} else {
		throw std::invalid_argument( "[error] Must specify plasma annulus with either PlasmaRadiusMin / PlasmaRadiusMax or AxialGapDistance and PlasmaColumnWidth" );
	}

	// Imposed Voltage
	readParameterFromFile(batch, "Voltage", ImposedVoltageVals, false, 0.0);

	// Wall Radius
	readParameterFromFile(batch, "WallRadius", WallRadiusVals, true);

	// Plasma Length
	readParameterFromFile(batch,"PlasmaLength", PlasmaLengthVals, true);

	// Auxiliary Heating
	readParameterFromFile(batch, "AuxiliaryHeating", AuxiliaryHeatingVals, false, 0.0);

	// This overrides the default for the chosen fuel if specified in file
	if ( batch.count( "IncludeAlphaHeating" ) == 1 )
	{
		bool aHeating = batch.at( "IncludeAlphaHeating" ).as_boolean();
		if ( aHeating )
			IncludeAlphaHeating = true;
		else
			IncludeAlphaHeating = false;
	}
	else
		IncludeAlphaHeating = std::nullopt;

	// This overrides the default for the chosen fuel
	if ( batch.count( "ReportNuclearDiagnostics" ) == 1 )
	{
		bool nucDiagnostics = batch.at( "ReportNuclearDiagnostics" ).as_boolean();
		if ( nucDiagnostics ) ReportNuclearDiagnostics = true;
		else ReportNuclearDiagnostics = false;
	}
	else ReportNuclearDiagnostics = std::nullopt;

	// Ion to electron temperature ratio
	readParameterFromFile(batch, "IonToElectronTemperatureRatio", TiTeVals, false, 0.0, true);

	// Effective Charge (Zeff)
	readParameterFromFile(batch, "Zeff", ZeffVals, false, 1.0, true);

	// Electron Density
	readParameterFromFile(batch, "ElectronDensity", ElectronDensityVals, true, 0.0, true);
	
	// Electron Temperature
	readParameterFromFile(batch, "ElectronTemperature", ElectronTemperatureVals, false, -1.0, false);

	// Neutral Density
	readParameterFromFile(batch, "NeutralDensity", NeutralDensityVals, false, 0.0, true);

	// Voltage Trace
	if ( batch.count( "VoltageTrace" ) == 1 ) {
		VoltageTrace = batch.at( "VoltageTrace" ).as_string();
		isTimeDependent = true;
	} else {
		VoltageTrace = "";
	}

	if ( batchConfig.count( "timestepping" ) == 1 ) {
		const auto & timestep_conf = toml::find<toml::table>( batchConfig, "timestepping" );
		OutputCadence = timestep_conf.at( "OutputCadence" ).as_floating();
		EndTime = timestep_conf.at( "EndTime" ).as_floating();
	} else {
		OutputCadence = 0.0005;
		EndTime = 1.00;
	}
}

void BatchRunner::cartesianProduct(vecMap& vectorOfMaps, mapSD& currentMap, vecPairVS::const_iterator currentI, vecPairVS::const_iterator end)
{
	if ( currentI == end )
	{
		vectorOfMaps.push_back(currentMap);
		return;
	}

	const std::pair<std::vector<double>*,std::string>& infoPair = *currentI;

	for ( std::vector<double>::const_iterator it = infoPair.first->begin() ; it != infoPair.first->end(); ++it )
	{
		currentMap.insert(std::pair<std::string, double>(infoPair.second,*it));
		cartesianProduct(vectorOfMaps,currentMap,currentI+1,end);
		currentMap.erase(infoPair.second);
	}
}

void BatchRunner::runBatchSolve()
{
	std::vector<std::map<std::string, double>> vectorOfMaps;
	std::map<std::string, double> currentMap;
	cartesianProduct(vectorOfMaps, currentMap, ptrsAndNamesToVectors.begin(), ptrsAndNamesToVectors.end());

	totalRuns = vectorOfMaps.size();

	if ( totalRuns > 1 && OutputFile == "" )
		throw std::invalid_argument("[error] Output file name is needed when running a batch solve");

	if ( totalRuns > 1 && isTimeDependent )
		throw std::invalid_argument( "[error] Multiple simultaneous time-dependent runs is not currently supported" );

	#pragma omp parallel for
	for ( int n = 0; n < totalRuns; n++ )
	{
		SolveIndividualMirrorPlasma(vectorOfMaps[n], n);
	}
#ifdef DEBUG
	std::cerr << "Total cases run:" << vectorOfMaps.size() << std::endl;
#endif
}

const double BatchRunner::step(std::vector<double> array)
{
	if ( array.size() != 3 || array[1] < array[0]) throw std::invalid_argument( "[error] Input must be configured [min, max, step size]" );
	if ( array[2] <= 0.0 ) return array[1] - array[0] + 1.0; // negative step implies only accpeting the first value of the matix
	else return array[2];
}

void BatchRunner::readParameterFromFile(toml::value batch, std::string configName, std::vector<double>& parameterVector, bool mandatory, double defaultValue, bool strictlyPositive )
{
	if ( batch.count( configName ) == 1 )
	{
		auto parameterData = toml::find(batch, configName );
		if ( parameterData.is_array() )
		{
			auto parameterArray = toml::get<std::vector<double>>(parameterData);
			for ( auto val = parameterArray[0]; val <= parameterArray[1]; val += step(parameterArray)) 
				parameterVector.push_back(val);
		}
		else if ( parameterData.is_floating() ) parameterVector.push_back(parameterData.as_floating());
		else if ( parameterData.is_integer() ) parameterVector.push_back(static_cast<double>(parameterData.as_floating()));
		else throw std::invalid_argument("[error] Non-compatible data given in configuration file");

		if ( strictlyPositive && *std::min_element( parameterVector.begin(), parameterVector.end() ) < 0.0 )
		throw std::invalid_argument( "[error] " + configName + " cannot be negative" );

	}
	else if ( !mandatory && batch.count( configName ) == 0 ) parameterVector.push_back(defaultValue);
	else if ( batch.count( configName ) > 1 ) throw std::invalid_argument("[error] " + configName + " cannot be specified more than once");
	else throw std::invalid_argument( "[error] " + configName + " unspecified or specified incorrectly");

	ptrsAndNamesToVectors.push_back(std::make_pair(&parameterVector,configName));
}

void BatchRunner::SolveIndividualMirrorPlasma(std::map<std::string, double> parameterMap, int currentRun)
{
	std::shared_ptr< MirrorPlasma > pReferencePlasmaState = std::make_shared<MirrorPlasma>( parameterMap,FuelName,ReportThrust,IncludeAlphaHeating,ReportNuclearDiagnostics, AmbipolarPhi, Collisional, IncludeCXLosses, OutputFile, NetcdfOutputFile, VoltageTrace);
	
	MCTransConfig config(pReferencePlasmaState, OutputCadence, EndTime);
	
	try
	{
		std::shared_ptr<MirrorPlasma> result = config.Solve();
		result->PrintReport(&parameterMap, currentRun, totalRuns);
		result->WriteNetCDFReport( &parameterMap, currentRun, totalRuns );
	}
	catch (int e ){}
	catch(boost::exception const&  ex){}
	catch(std::exception e){}
}

template<typename K, typename V>
void BatchRunner::print_maps(std::vector<std::map<K, V>> const &vec)
{
	for ( auto const &m: vec )
	{
		for ( auto const &pair: m ) 
			std::cout << "{" << pair.first << ": " << pair.second << "}\n";
		std::cout << std::endl;
	}
}
