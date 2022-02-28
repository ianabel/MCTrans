#include <toml.hpp>
#include <iostream>

#include "BatchRunner.hpp"
#include "MirrorPlasma.hpp"
#include "Config.hpp"

BatchRunner::BatchRunner(std::string const& batchFile)
{
    const auto batchConfig = toml::parse( batchFile );
    //Algorithm Parameters
    if(batchConfig.count( "algorithm" ) != 0)
    {
        const auto algConfig = toml::find<toml::table>( batchConfig, "algorithm" );

        readParameterFromFile(algConfig, "ParallelFudgeFactor", ParallelFudgeFactorVals, false, 1.0);
        readParameterFromFile(algConfig, "PerpFudgeFactor", PerpFudgeFactorVals, false, 1.0);
        readParameterFromFile(algConfig, "InitialTemp", InitialTempVals, false, 0.1);
        readParameterFromFile(algConfig, "InitialMach", InitialMachVals, false, 4.0);
        
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
		Collisional = false;
    }
    else
    {
        ParallelFudgeFactorVals.push_back(1.0);
		PerpFudgeFactorVals.push_back(1.0);
		InitialTempVals.push_back(0.1);
		InitialMachVals.push_back(4.0);
		AmbipolarPhi = true;
		Collisional = false;
		OutputFile  = "";
		NetcdfOutputFile = "";
    }

    const auto batch = toml::find<toml::value>( batchConfig, "configuration" );
    //Fuel name
    if ( batch.count( "IonSpecies" ) != 1 )
		throw std::invalid_argument( "Fuel must be specified once in the [configuration] block" );
	FuelName = batch.at( "IonSpecies" ).as_string();

    //Report Thrust
    if ( batch.count( "ReportThrust" ) == 1 )
		ReportThrust = batch.at( "ReportThrust" ).as_boolean();
	else
		ReportThrust = false;

    //Central Field Strength
    readParameterFromFile(batch, "CentralCellField", CentralCellFieldStrengthVals);
    
    //Mirror Ratio
    if ( batch.count( "MirrorRatio" ) == 1 ) 
    {
        readParameterFromFile(batch, "MirrorRatio", MirrorRatioVals, false);
        useMirrorRatio = true;
        if ( batch.count( "ThroatField" ) == 1 ) {
			throw std::invalid_argument( toml::format_error( "[error] Cannot specify botth Miror Ratio and Throat Field",batch.at( "MirrorRatio" )," mirror ratio defined here", batch.at( "ThroatField" ), " Throat field here" ) );
		}
    }
    //Magnetic Throat Field
	else if ( batch.count( "ThroatField" ) == 1 ) 
    {
		readParameterFromFile(batch, "ThroatField", MagFieldThroatVals, false);
        useMirrorRatio = false;
		if ( batch.count( "MirrorRatio" ) == 1 ) {
			throw std::invalid_argument( toml::format_error( "[error] Cannot specify botth Miror Ratio and Throat Field",batch.at( "MirrorRatio" )," mirror ratio defined here", batch.at( "ThroatField" ), " Throat field here" ) );
		}
	} 
    else {
		throw std::invalid_argument( "[error] Must specify either MirrorRatio or ThroatField" );
	}

    //Plasma Radius Max and Min
    if ( batch.count( "PlasmaRadiusMin" ) == 1  || batch.count( "PlasmaRadiusMax" ) == 1 ) {
		if ( batch.count( "PlasmaRadiusMin" ) == 0 )
			throw std::invalid_argument( toml::format_error( "[error] When PlasmaRadiusMax is specified you must also set PlasmaRadiusMin",batch.at( "PlasmaRadiusMax" )," max radius set here" ) );
		if ( batch.count( "PlasmaRadiusMax" ) == 0 )
			throw std::invalid_argument( toml::format_error( "[error] When PlasmaRadiusMin is specified you must also set PlasmaRadiusMax",batch.at( "PlasmaRadiusMin" )," min radius set here" ) );
        
        readParameterFromFile(batch, "PlasmaMinRadius", AxialGapDistanceVals, false);
        readParameterFromFile(batch, "PlasmaMaxRadius", PlasmaMaxRadiusVals, false);
        usePlasmaRadiusMaxMin = true;

		if ( batch.count( "AxialGapDistance" ) == 1 ) {
			throw std::invalid_argument( toml::format_error( "[error] When PlasmaRadiusMin / Max are specified you cannot set AxialGapDistance",batch.at( "PlasmaRadiusMin" )," minimum radius set here", batch.at( "AxialGapDistance" ), " AxialGapDistance found here" ) );
		} else if ( batch.count( "PlasmaColumnWidth") == 1 ) {
			throw std::invalid_argument( toml::format_error( "[error] When PlasmaRadiusMin / Max are specified you cannot set AxialGapDistance",batch.at( "PlasmaRadiusMin" )," minimum radius set here", batch.at( "AxialGapDistance" ), " AxialGapDistance found here" ) );
		}
    }

    //Axial Gap Distance and Plasma Column Width
	else if ( batch.count( "AxialGapDistance" ) == 1  || batch.count( "PlasmaColumnWidth" ) == 1 ) {
		if ( batch.count( "AxialGapDistance" ) == 0 )
			throw std::invalid_argument( toml::format_error( "[error] When PlasmaColumnWidth is specified you must also set AxialGapDistance",batch.at( "PlasmaColumnWidth" )," width set here" ) );
		if ( batch.count( "PlasmaColumnWidth" ) == 0 )
			throw std::invalid_argument( toml::format_error( "[error] When AxialGapDistance is specified you must also set PlasmaColumnWidth",batch.at( "AxialGapDistance" )," axial gap set here" ) );
		
        readParameterFromFile(batch, "AxialGapDistance", AxialGapDistanceVals, false);
        readParameterFromFile(batch,"PlasmaColumnWidth", PlasmaColumnWidthVals, false);
        usePlasmaRadiusMaxMin = false;

		if ( batch.count( "PlasmaMinRadius" ) == 1 ) {
			throw std::invalid_argument( toml::format_error( "[error] When AxialGapDistance and PlasmaColumnWidth are specified you cannot set PlasmaMinRadius",batch.at( "PlasmaRadiusMin" )," minimum radius set here", batch.at( "AxialGapDistance" ), " AxialGapDistance found here" ) );
		} else if ( batch.count( "PlasmaMaxRadius") == 1 ) {
			throw std::invalid_argument( toml::format_error( "[error] When AxialGapDistance and PlasmaColumnWidth are specified you cannot set PlasmaMaxRadius",batch.at( "PlasmaRadiusMax" )," maximum radius set here", batch.at( "AxialGapDistance" ), " AxialGapDistance found here" ) );
		}
	} else {
		throw std::invalid_argument( "[error] Must specify plasma annulus with either PlasmaRadiusMin / PlasmaRadiusMax or AxialGapDistance and PlasmaColumnWidth" );
	}

    //Imposed Voltage
    readParameterFromFile(batch, "Voltage", ImposedVoltageVals, false, 0.0);

    //Wall Radius
    readParameterFromFile(batch, "WallRadius", WallRadiusVals, true);

    //Plasma Length
    readParameterFromFile(batch,"PlasmaLength", PlasmaLengthVals, true);

    //Auxiliary Heating
    readParameterFromFile(batch, "AuxiliaryHeating", AuxiliaryHeatingVals, false, 0.0);

    // This overrides the default for the chosen fuel if specified in file
	if ( batch.count( "IncludeAlphaHeating" ) == 1 )
    {
		bool aHeating = batch.at( "IncludeAlphaHeating" ).as_boolean();
        if(aHeating) AlphaHeating = tru;
        else AlphaHeating = fal;
    }
    else AlphaHeating = unspecified;

	// This overrides the default for the chosen fuel
	if ( batch.count( "ReportNuclearDiagnostics" ) == 1 )
    {
		bool nucDiagnostics = batch.at( "ReportNuclearDiagnostics" ).as_boolean();
        if(nucDiagnostics) ReportNuclearDiagnostics = tru;
        else ReportNuclearDiagnostics = fal;
    }
    else ReportNuclearDiagnostics = unspecified;

    //Ion to electron temperature ratio
    readParameterFromFile(batch, "IonToElectronTemperatureRatio", TiTeVals, false, 0.0, true);

    //Effective Charge (Zeff)
    readParameterFromFile(batch, "Zeff", ZeffVals, false, 1.0, true);

    //Electron Density
    readParameterFromFile(batch, "ElectronDensity", ElectronDensityVals, true, 0.0, true);
    
    //Electron Temperature
    readParameterFromFile(batch, "ElectronTemperature", ElectronTemperatureVals, false, -1.0, false);

    //Neutral Density
    readParameterFromFile(batch, "NeutralDensity", NeutralDensityVals, false, 0.0, false);

    //Voltage Trace
    if ( batch.count( "VoltageTrace" ) == 1 ) {
		VoltageTrace = batch.at( "VoltageTrace" ).as_string();
	} else VoltageTrace = "";
}

void BatchRunner::cartesianProduct(std::vector<std::map<std::string, double>>& vectorOfMaps, std::map<std::string, double>& currentMap, std::vector<std::pair<std::vector<double>*,std::string>>::const_iterator currentI, std::vector<std::pair<std::vector<double>*,std::string>>::const_iterator end)
{
    if(currentI == end)
    {
        vectorOfMaps.push_back(currentMap);
        return;
    }

    const std::pair<std::vector<double>*,std::string>& infoPair = *currentI;

    for(std::vector<double>::const_iterator it = infoPair.first->begin() ; it != infoPair.first->end(); ++it)
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

    //print_maps(vectorOfMaps);

    for(auto caseMap : vectorOfMaps)
    {
        SolveIndividualMirrorPlasma(caseMap);
    }

    std::cout << "Total cases run:" << vectorOfMaps.size() << std::endl;
}

const double BatchRunner::step(std::vector<double> array)
{
    if(array.size() != 3 || array[1] < array[0]) throw std::invalid_argument( "Input must be configured [min, max, step size]" );
    if(array[2] <= 0.0) return array[1] - array[0] + 1.0; // negative step implies only accpeting the first value of the matix
    else return array[2];
}

void BatchRunner::readParameterFromFile(toml::value batch, std::string configName, std::vector<double>& parameterVector, bool mandatory, double defaultValue, bool strictlyPositive )
{
    if( batch.count( configName ) == 1)
    {
        auto parameterData = toml::find(batch, configName );
        if(parameterData.is_array())
        {
            auto parameterArray = toml::get<std::vector<double>>(parameterData);
            for(auto val = parameterArray[0]; val <= parameterArray[1]; val += step(parameterArray)) 
                parameterVector.push_back(val);
        }
        else if(parameterData.is_floating()) parameterVector.push_back(parameterData.as_floating());
        else if(parameterData.is_integer()) parameterVector.push_back(static_cast<double>(parameterData.as_floating()));
        else throw std::invalid_argument("Non-compatible data given in configuration file");

        if(strictlyPositive && *std::min_element( parameterVector.begin(), parameterVector.end() ) < 0.0 )
        throw std::invalid_argument( configName + "cannot be negative" );

    }
    else if(!mandatory && batch.count( configName ) == 0) parameterVector.push_back(defaultValue);
    else if( batch.count( configName ) > 1 ) throw( configName + "cannot be specified more than once");
    else throw( configName + "unspecified or specified incorrectly");

    ptrsAndNamesToVectors.push_back(std::make_pair(&parameterVector,configName));
}

void BatchRunner::SolveIndividualMirrorPlasma(std::map<std::string, double> parameterMap)
{
    std::shared_ptr< MirrorPlasma::VacuumMirrorConfiguration > pVacuumConfig = std::make_shared<MirrorPlasma::VacuumMirrorConfiguration>( parameterMap,FuelName,reportThrust,AlphaHeating,ReportNuclearDiagnostics, AmbipolarPhi, Collisional, OutputFile, NetcdfOutputFile );
    std::shared_ptr< MirrorPlasma > pReferencePlasmaState = std::make_shared<MirrorPlasma>(pVacuumConfig, parameterMap, VoltageTrace);
    MCTransConfig config(pReferencePlasmaState);

    std::shared_ptr<MirrorPlasma> result = config.Solve();
	result->PrintReport();
}

int BatchRunner::totalCases()
{
    int total = CentralCellFieldStrengthVals.size()*
    AxialGapDistanceVals.size()*
    ImposedVoltageVals.size()*
    PlasmaLengthVals.size()*
    WallRadiusVals.size()*
    AuxiliaryHeatingVals.size()*
    ZeffVals.size()*
    ElectronDensityVals.size()*
    ElectronTemperatureVals.size()*
    TiTeVals.size();

    if(useMirrorRatio) total *= MirrorRatioVals.size();
    else total *= MagFieldThroatVals.size();

    if(usePlasmaRadiusMaxMin) total *= PlasmaMaxRadiusVals.size();
    else total *= PlasmaColumnWidthVals.size();

    //Handy debugging comment block
    /*
    std::cout << "Central Field Strength    " << CentralCellFieldStrengthVals.size() << std::endl;
    std::cout << "Throat Field Strength    " << MagFieldThroatVals.size() << std::endl;
    std::cout << "Mirror Ratio    " << MirrorRatioVals.size() << std::endl;
    std::cout << "Plasma Max Radius    " << PlasmaMaxRadiusVals.size() << std::endl;
    std::cout << "Axial Gap Distance    " << AxialGapDistanceVals.size() << std::endl;
    std::cout << "Plasma Column Width    " << PlasmaColumnWidthVals.size() << std::endl;
    std::cout << "Imposed Voltage    " << ImposedVoltageVals.size() << std::endl;
    std::cout << "Plasma Length    " << PlasmaLengthVals.size() << std::endl;
    std::cout << "Wall Radius    " << WallRadiusVals.size() << std::endl;
    std::cout << "Auxiliary Heating    " << AuxiliaryHeatingVals.size() << std::endl;
    std::cout << "Zeff    " << ZeffVals.size() << std::endl;
    std::cout << "Electron Density    " << ElectronDensityVals.size() << std::endl;
    std::cout << "Electron Temperature    " << ElectronTemperatureVals.size() << std::endl;
    std::cout << "TiTe    " << TiTeVals.size() << std::endl;
    */

    return total;
}

template<typename K, typename V>
void BatchRunner::print_maps(std::vector<std::map<K, V>> const &vec)
{
    for(auto const &m: vec)
    {
        for (auto const &pair: m) 
            std::cout << "{" << pair.first << ": " << pair.second << "}\n";
        std::cout << std::endl;
    }
}