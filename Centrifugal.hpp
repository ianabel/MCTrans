#ifndef CENTRIFUGAL_HPP
#define CENTRIFUGAL_HPP

#include <string>

struct Plasma {
	double Mu; // Units of the proton mass
	double Zeff; // For Bremsstrahlung
	double ZIon; // For cyclotron frequencies etc.
	double TiTe; // T_i / T_e
	double ElectronTemperature; // Units of keV
	double ElectronDensity; // Units of 10^20
	std::string Fuel;
	bool LumpedImpurity; // Explicit impurity species?
	double IonDensity;
	double ImpurityDensity;
	double NeutralDensity;
	double NeutralSource;
};

struct Configuration {
	double MirrorRatio;
	double PlasmaColumnWidth;
	double PlasmaLength;
	double AxialGapDistance;
	double WallRadius;
	double CentralCellFieldStrength;
	double MachNumber;
	double Voltage;
	double AuxiliaryHeating; // in MW
	bool ReportNuclear;
	bool IncludeAlphaHeating;
	bool MachSolve;
	bool TempSolve;
};

// SI Units, except with k_B = 1 and temperature in energy
// units
constexpr double ElectronMass = 9.10938356e-31; 
constexpr double VacuumPermittivity = 8.8541878128e-12;
constexpr double PermeabilityOfFreeSpace = 1.25663706212e-6; // post-2018 not exactly 4*pi * 10^-7
constexpr double ElectronCharge = 1.60217662e-19;
constexpr double ProtonMass = 1.6726219e-27;

// Reference temperature of 1 keV, expressed in Joules
constexpr double ReferenceTemperature = 1000 * ElectronCharge;
constexpr double ReferenceDensity = 1e20;

// Reference Cross-section for Charge-Exchange neutral cross-section
// in cm^2
constexpr double CXCrossSection = 1e-15;
constexpr double IonizationCrossSection = 1e-16;


// Finds the steady-state Mach Number
void SetMachNumber( Plasma const&, Configuration & );
// Find T given M (or V)
void SetTemperature( Plasma &, Configuration const& );

double PlasmaVolume( Configuration const& );
double WallArea( Configuration const& );
double BremsstrahlungLosses( Plasma const& );
double AlfvenMachNumber( Plasma const&, Configuration const& );
double ViscousHeating( Plasma const& , Configuration const& );
double TauEInverse( Plasma const& , Configuration const& );
double ParallelElectronHeatLoss( Plasma const& , Configuration const& );
double ClassicalHeatLossRate( Plasma const& , Configuration const& );
double IonCollisionTime( Plasma const& );
double ElectricPotential( Plasma const &plasma, Configuration const& conf );
double IonElectronEquilibrationTime( Plasma const& );
double IonLarmorRadius( Plasma const& plasma, Configuration const& conf );
double Beta( Plasma const& plasma, Configuration const& conf );
double DebyeLength( Plasma const& plasma );
double CentrifugalPotential( Plasma const& plasma, Configuration const& conf );
double NuStar( Plasma const& plasma, Configuration const& conf );
double ParallelKineticEnergyLoss( Plasma const& plasma, Configuration const& conf );
double ThermalEnergy( Plasma const& plasma, Configuration const& conf );
double KineticEnergy( Plasma const& plasma, Configuration const& conf );


#endif // CENTRIFUGAL_HPP
