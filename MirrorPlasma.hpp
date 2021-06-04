#ifndef MIRRORPLASMA_HPP
#define MIRRORPLASMA_HPP

#include "PlasmaPhysics.hpp"
#include "Species.hpp"
#include "toml11/toml.hpp"

#include <memory>
#include <cmath>

class MirrorPlasma {
	public:
		class VacuumMirrorConfiguration {
			public:
				VacuumMirrorConfiguration( toml::value const& );
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

				bool AlphaHeating;
				bool ReportNuclearDiagnostics;
				bool ReportMomentumLoss;

				double PlasmaVolume() const {
					return M_PI * ( PlasmaColumnWidth + 2 * AxialGapDistance ) * PlasmaColumnWidth * PlasmaLength;
				};
				double WallArea() const {
					return 2.0 * M_PI * WallRadius * WallRadius * PlasmaLength; 
				};

				double ImposedVoltage;
			private:

		};

		// Copy Constructor
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

		MirrorPlasma( toml::value const& configSection );

		double FuellingRate;
		double IonDensity,IonTemperature;
		double ElectronDensity,ElectronTemperature;
		double NeutralDensity,NeutralSource;
		double Zeff;

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

		void PrintReport() const;

		std::shared_ptr< VacuumMirrorConfiguration > pVacuumConfig;

		void SetMachFromVoltage();
		double AmbipolarPhi() const;

		void ComputeSteadyStateNeutrals();

		double ParallelMomentumLossRate() const;
	private:

		double LogLambdaElectron() const;
		double LogLambdaIon() const;
		double LogLambdaAlphaElectron() const;


		double ElectronCollisionTime() const;
		double IonCollisionTime() const;
		double CollisionalTemperatureEquilibrationTime() const;

		double SlowingDownTime() const;

		double IonToElectronHeatTransfer() const;

		double BremsstrahlungLosses() const;
		double NeutralLosses() const;
		double ClassicalHeatLosses() const;
		double ParallelHeatLosses() const;

		double ParallelElectronPastukhovLossRate( double Phi ) const;
		double ParallelElectronParticleLoss() const;
		double ParallelElectronHeatLoss() const;

		double ParallelIonPastukhovLossRate( double Phi ) const;
		double ParallelIonParticleLoss() const;
		double ParallelIonHeatLoss() const;

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
		double NeutronOutput() const;
		double ThermalPowerOutput() const;
		double NeutronWallLoading() const;
		double DDNeutronRate() const;

		double AlphaHeating() const;
		double PromptAlphaLossFraction() const;
		double AlphaPromptLosses() const;
		double AlphaParallelLossRate() const;
};




#endif // MIRRORPLASMA_HPP

