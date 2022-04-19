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
	double Mass; // kg
	std::string Name; // For reporting
} Species;

#endif // Definition of species struct
