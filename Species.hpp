#ifndef SPECIES_HPP
#define SPECIES_HPP

#include <string>
#include "PlasmaPhysics.hpp"

typedef struct Species_t {
	enum Type {
		Electron,
		Ion,
		TraceImpurity,
		Neutral
	} type;
	double Charge; // Units of e
	double Mass; // Units of Proton Mass
	std::string Name; // For reporting
} Species;

extern Species Electron;
extern Species Proton;
extern Species Deuteron;
extern Species NeutralHydrogen;


#endif // Definition of species struct
