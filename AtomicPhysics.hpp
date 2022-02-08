
#include "Species.hpp"
#include <functional>

class CrossSection {
	public:
		CrossSection( std::function<double( double )> sigma, double E_min, double E_max, Species s1, Species s2 )
		{
			sigmaImplementation = sigma;
			MinEnergy = E_min;
			MaxEnergy = E_max;
			Species1 = s1;
			Species2 = s2;
			ReducedMass = Species1.Mass * Species2.Mass / ( Species1.Mass + Species2.Mass );
		};

		// Masses in units of the proton mass
		double ReducedMass;
		// in units of electron volts (not keV!)
		double MinEnergy,MaxEnergy;
		// Species types
		Species_t Species1,Species2;

		// Return the cross section in cm^2 -- not barns, not m^2
		// takes Centre of Mass energy in electron Volts.
		double operator()( double CentreOfMassEnergy ) const
		{
			// Should we check E_min/E_max bounds here?
			return sigmaImplementation( CentreOfMassEnergy );
		}

	private:
		std::function<double( double )> sigmaImplementation;
};

// From Janev (1993) & Janev (1987), argument is centre-of-mass energy in electron Volts, 
// return values are in cm^2
double protonImpactIonizationCrossSection( double );
double HydrogenChargeExchangeCrossSection( double );
double electronImpactIonizationCrossSection( double );

// Cross Section objects for integration
CrossSection protonImpactIonization( protonImpactIonizationCrossSection, 200, 1e6, Proton, NeutralHydrogen );
CrossSection HydrogenChargeExchange( HydrogenChargeExchangeCrossSection, 0.1, 1e6, Proton, NeutralHydrogen );
CrossSection electronImpactIonization( electronImpactIonizationCrossSection, 13.6, 1e6, Electron, NeutralHydrogen );

double neutralsRateCoefficentHot( CrossSection const & , std::shared_ptr<const MirrorPlasma> );
double neutralsRateCoefficentCold( CrossSection const &, std::shared_ptr<const MirrorPlasma> );

