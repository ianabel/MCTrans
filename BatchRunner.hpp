#pragma once
// This class is used to read batch cases into the solver and run many cases. Each value can be assigned
// a 3 value array [min, max, stepsize]. This class creates vectors containing all the values to be run for
// each parameter. n^m mirrorPlasmas are created and each is run as a seperate case.
// Note you can also specify a single value to be read and this just gets taken as a single input for that variable
#include <string>
#include <vector>
#include <toml.hpp>
#include <optional>

#include "MirrorPlasma.hpp"

class BatchRunner {
public:
	BatchRunner(std::string const& batchFile);
	void runBatchSolve();
private:
	// solving function. Builds the mirrorplasma object, excecutes a solve and prints the report
	void SolveIndividualMirrorPlasma(std::map<std::string, double> parameterMap, int currentRun);

	// Catchall function which takes parameters from the toml object and populates the parameter vectors
	void readParameterFromFile(toml::value batch, std::string configName, std::vector<double>& parameterVector, bool mandatory = true, double defaultValue = 0.0, bool strictlyPositive = false);
	
	// calculates the necessary increment in a parameter at each step
	const double step(std::vector<double> array);

	// recursive function which takes creates a map of each mirror plasma case and stores all the maps in a vector
	void cartesianProduct(std::vector<std::map<std::string, double>>& final, std::map<std::string, double>& current, std::vector<std::pair<std::vector<double>*,std::string>>::const_iterator currentInput, std::vector<std::pair<std::vector<double>*,std::string>>::const_iterator end);
	std::string FuelName;

	template<typename K, typename V>
	void print_maps(std::vector<std::map<K, V>> const &vec);
	
	std::vector<std::pair<std::vector<double>*,std::string>> ptrsAndNamesToVectors;

	// VacuumConfig parameters
	std::vector<double> CentralCellFieldStrengthVals;
	std::vector<double> MagFieldThroatVals;
	std::vector<double> MirrorRatioVals;
	std::vector<double> PlasmaMaxRadiusVals;
	std::vector<double> AxialGapDistanceVals;
	std::vector<double> PlasmaColumnWidthVals;
	std::vector<double> ImposedVoltageVals;
	std::vector<double> PlasmaLengthVals;
	std::vector<double> WallRadiusVals;
	std::vector<double> AuxiliaryHeatingVals;
	std::vector<double> NeutralDensityVals;

	std::optional<bool> IncludeAlphaHeating = std::nullopt;
	std::optional<bool> ReportNuclearDiagnostics = std::nullopt;
	bool ReportThrust;

	// MirrorPlasma parameters
	std::vector<double> ZeffVals;
	std::vector<double> ElectronDensityVals,ElectronTemperatureVals;
	std::vector<double> TiTeVals;

	// Algorithm parameters
	std::vector<double> ParallelFudgeFactorVals;
	std::vector<double> PerpFudgeFactorVals;
	std::vector<double> InitialTempVals;
	std::vector<double> InitialMachVals;
	std::vector<double> RateThresholdVals;
	std::vector<double> SundialAbsTolVals;
	std::vector<double> SundialsRelTolVals;

	bool AmbipolarPhi = true;
	bool Collisional = false;
	bool IncludeCXLosses = false;
	std::string OutputFile  = "";
	std::string NetcdfOutputFile = "";

	// Time dependent parameters
	std::string VoltageTrace;
	double EndTime,OutputCadence;
	bool isTimeDependent = false;

	int totalRuns = 1;
};
