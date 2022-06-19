#ifndef MIRRORPLASMA_HPP
#define MIRRORPLASMA_HPP

#include "PlasmaPhysics.hpp"
#include "Species.hpp"
#include "NetCDFIO.hpp"
#include <toml.hpp>

#include <boost/math/interpolators/barycentric_rational.hpp>

#include <memory>
#include <cmath>
#include <map>

class MirrorPlasma {
	public:
		class VacuumMirrorConfiguration {
			public:
				// VacuumMirrorConfiguration( toml::value const& );
				VacuumMirrorConfiguration(const std::map<std::string, double>& parameterMap, std::string FuelName, bool reportThrust, std::optional<bool> AlphaHeating, std::optional<bool> ReportNuclearDiagnostics, bool ambiPolPhi, bool collisions, bool includeCXLosses, std::string asciiOut, std::string netCdfOut);
				/*
				VacuumMirrorConfiguration( const& VacuumMirrorConfiguration other ) :
					IonSpecies( other.IonSpecies )
				{
					MirrorRatio = other.MirrorRatio;
					PlasmaColumnWidth = other.PlasmaColumnWidth;
					PlasmaLength = other.PlasmaLength;
					AxialGapDistance = other.AxialGapDistance;
					WallRadius = other.WallRadius;
					CentralCellFieldStrength = other.CentralCellFieldStrength;
					AuxiliaryHeating = other.AuxiliaryHeating;
				}
				*/
				Species_t IonSpecies;
				double MirrorRatio;
				double PlasmaColumnWidth;
				double PlasmaLength;
				double AxialGapDistance;
				double WallRadius;
				double CentralCellFieldStrength;
				double AuxiliaryHeating;

				double ParallelFudgeFactor;
				double PerpFudgeFactor;
				bool AmbipolarPhi;

				bool IncludeCXLosses;

				bool AlphaHeating;
				bool ReportNuclearDiagnostics;
				bool ReportThrust;

				double InitialTemp;
				double InitialMach;
				bool Collisional;

				std::string OutputFile;
				std::string NetcdfOutputFile;

				double SundialsAbsTol,SundialsRelTol;
				double RateThreshold;

				double PlasmaVolume() const {
					return M_PI * ( PlasmaColumnWidth + 2 * AxialGapDistance ) * PlasmaColumnWidth * PlasmaLength;
				};
				double WallArea() const {
					return 2.0 * M_PI * WallRadius * WallRadius * PlasmaLength;
				};

				double ImposedVoltage;
				double PlasmaInnerRadius() const { return AxialGapDistance; };
				double PlasmaOuterRadius() const { return AxialGapDistance + PlasmaColumnWidth; };
				double PlasmaCentralRadius() const { return AxialGapDistance + PlasmaColumnWidth / 2.0; };
			private:

		};

		// Copy Constructor
		/*
		MirrorPlasma( MirrorPlasma const& other ) :
			pVacuumConfig( other.pVacuumConfig )
		{
			FuellingRate = other.FuellingRate;
			IonDensity = other.IonDensity;
			IonTemperature = other.IonTemperature;
			ElectronDensity = other.ElectronDensity;
			ElectronTemperature = other.ElectronTemperature;
		   NeutralDensity = other.NeutralDensity;
			NeutralSource = other.NeutralSource;
			Zeff = other.Zeff;
		   MachNumber = other.MachNumber;
		};
		*/

		MirrorPlasma( MirrorPlasma const& ) = delete;

		// MirrorPlasma( toml::value const& configSection );

		MirrorPlasma(std::shared_ptr< VacuumMirrorConfiguration > pVacuumConfig, std::map<std::string,double> parameterMap, std::string vTrace);

		double FuellingRate;
		double IonDensity,IonTemperature;
		double ElectronDensity,ElectronTemperature;
		double NeutralDensity,NeutralSource;
		double Zeff;

		bool FixedNeutralDensity;

		double MachNumber; // Sonic Mach number defined with c_s^2 = Z T_e/m_i

		double SoundSpeed() const {
			// We *define* c_s^2 = Z_i T_e / m_i.
			double cs = ::sqrt( pVacuumConfig->IonSpecies.Charge * ElectronTemperature*ReferenceTemperature / ( pVacuumConfig->IonSpecies.Mass * ProtonMass ) );
			return cs;
		};

		// Defining (3/2) d( n_i T_i )/dt = IonHeating - IonHeatLosses
		//
		double ElectronHeatLosses() const;
		double IonHeatLosses() const;
		double ElectronHeating() const;
		double IonHeating() const;


		double EnergyConfinementTime() const {
			double StoredEnergy = 1.5 * ( ElectronDensity * ElectronTemperature + IonDensity * IonTemperature ) * ReferenceDensity * ReferenceTemperature;
			double TotalHeatLosses = ElectronHeatLosses() + IonHeatLosses();
			return StoredEnergy / TotalHeatLosses;
		};

		// RHS of dn/dt = <stuff>
		double IonParticleLosses() const;
		double ElectronParticleLosses() const;

		double ElectricPotential() const {
			return pVacuumConfig->CentralCellFieldStrength * MachNumber * SoundSpeed() * pVacuumConfig->PlasmaColumnWidth;
		};
		double IonLarmorRadius() const {
			return 1.02 * ::sqrt( pVacuumConfig->IonSpecies.Mass * IonTemperature * 1000 ) / ( pVacuumConfig->IonSpecies.Charge * pVacuumConfig->CentralCellFieldStrength * 10000 );
		};

		double Beta() const;
		double DebyeLength() const;

		double NuStar() const;
		double KineticEnergy() const;
		double ThermalEnergy() const;

		void PrintReport(std::map<std::string, double>* parameterMap = nullptr, int currentRun = 1, int totalRun = 1);
		void WriteNetCDFReport(std::map<std::string, double>* parameterMap = nullptr, int currentRun = 1, int totalRun = 1);

		void InitialiseNetCDF();
		void WriteTimeslice( double T );
		void FinaliseNetCDF();


		std::shared_ptr< VacuumMirrorConfiguration > pVacuumConfig;

		void SetMachFromVoltage();
		double AmbipolarPhi() const;

		void ComputeSteadyStateNeutrals();

		double TotalAngularMomentumLosses() const;
		double MomentOfInertia() const;
		double InjectedTorque( double ) const;
		double ParallelMomentumLossRate() const;
		double initialTemperature() const { return pVacuumConfig->InitialTemp; };
		double initialMach() const { return pVacuumConfig->InitialMach; };

		bool isSteady;
	private:
		NetCDFIO nc_output;

		double LogLambdaElectron() const;
		double LogLambdaIon() const;
		double LogLambdaAlphaElectron() const;


		double ElectronCollisionTime() const;
		double IonCollisionTime() const;
		double CollisionalTemperatureEquilibrationTime() const;

		double SlowingDownTime() const;

		double IonToElectronHeatTransfer() const;

		double BremsstrahlungLosses() const;
		double CyclotronLosses() const;
		double RadiationLosses() const;

		double CXLossRate() const;
		double CXHeatLosses() const;
		double CXMomentumLosses() const;

		double ClassicalHeatLosses() const;
		double ParallelHeatLosses() const;

		double ParallelElectronPastukhovLossRate( double Phi ) const;
		double ParallelElectronParticleLoss() const;
		double ParallelElectronHeatLoss() const;

		double Chi_i( double ) const;
		double Chi_i() const;

		double ParallelIonPastukhovLossRate( double Phi ) const;
		double ParallelIonParticleLoss() const;
		double ParallelIonHeatLoss() const;

		double AngularMomentumPerParticle() const;

		double ParallelKineticEnergyLoss() const;

		double ClassicalIonHeatLoss() const;
		double ClassicalElectronHeatLoss() const;

		double ClassicalElectronParticleLosses() const;
		double ClassicalIonParticleLosses() const;

		double ClassicalViscosity() const;
		double AlfvenMachNumber() const;

		double ViscousHeating() const;

		double CentrifugalPotential() const;

		double IonCyclotronFrequency() const
		{
			double MagneticField = pVacuumConfig->CentralCellFieldStrength;
			return pVacuumConfig->IonSpecies.Charge * ElectronCharge * MagneticField / ( pVacuumConfig->IonSpecies.Mass * ProtonMass );
		};

		double ElectronCyclotronFrequency() const
		{
			double MagneticField = pVacuumConfig->CentralCellFieldStrength;
			return ElectronCharge * MagneticField / ElectronMass;
		};

		double FusionAlphaPowerDensity() const;
		double AlphaProductionRate() const;
		double NeutronOutput() const;
		double ThermalPowerOutput() const;
		double NeutronWallLoading() const;
		double DDNeutronRate() const;

		using interpolant = boost::math::barycentric_rational<double>;
		std::unique_ptr<interpolant> VoltageFunction;
		void ReadVoltageFile( std::string const& );
		double time;

	public:
		double AlphaHeating() const;
		double PromptAlphaLossFraction() const;
		double PromptAlphaThrust() const;
		double AlphaPromptLosses() const;
		double AlphaParallelLossRate() const;
		double ViscousTorque() const;
		double ParallelAngularMomentumLossRate() const;
		double RadialCurrent() const;
		double ParallelIonThrust() const;
		double ParallelCurrent(double) const;
		void UpdateVoltage();
		void SetTime( double );
		bool isTimeDependent;

		double ExternalResistance;
};




#endif // MIRRORPLASMA_HPP
