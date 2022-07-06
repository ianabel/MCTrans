#ifndef ATOMICPHYSICS_HPP
#define ATOMICPHYSICS_HPP

#include "Species.hpp"
#include <memory>
#include <functional>


class CrossSection {
	public:
		CrossSection( std::function<double( double )> sigma, double E_min, double E_max, Species particle, Species target )
		{
			sigmaImplementation = sigma;
			MinEnergy = E_min;
			MaxEnergy = E_max;
			Particle = particle;
			Target = target;
			ReducedMass = particle.Mass * target.Mass / ( particle.Mass + target.Mass );
		};

		// Masses in units of the proton mass
		double ReducedMass;
		// in units of electron volts (not keV!)
		double MinEnergy,MaxEnergy;
		// Species types
		Species Particle,Target;

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

double neutralsRateCoefficientHot( CrossSection const & , std::shared_ptr<const MirrorPlasma> );
double neutralsRateCoefficientCold( CrossSection const &, std::shared_ptr<const MirrorPlasma> );
#endif // ATOMICPHYSICS_HPP
