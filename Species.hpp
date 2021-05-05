#ifndef SPECIES_HPP
#define SPECIES_HPP

#include <string>

struct Species_t {
	enum Type {
		Electron,
		Ion,
		TraceImpurity
	} type;
	double Charge; // Units of e
	double Mass; // Units of Proton Mass
	std::string Name; // For reporting
};

#endif // Definition of species struct
