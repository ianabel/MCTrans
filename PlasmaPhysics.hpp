#ifndef PLASMAPHYSICS_HPP
#define PLASMAPHYSICS_HPP

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

// e in Gaussian Units is e_SI / sqrt( 4 * pi * epsilon_0 )
// precomuting this to e_SI / (sqrt(4*pi))*(sqrt(epsilon_0)) gives:
constexpr double GaussianElectronCharge = ElectronCharge / ( 3.54490770181103205458 * .00000297559873134802 );

// CODATA value 6.644657230(82)×10−27 kg
constexpr double AlphaMass = 6.644657230e-27;
constexpr double Z_Alpha = 2.0;

enum Fuel {
	Hydrogen,
	Deuterium,
	DTFuel,
};

#endif // PLASMAPHYSICS_HPP
