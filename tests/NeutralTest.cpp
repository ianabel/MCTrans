
#include "AtomicPhysics.hpp"
#include "MirrorPlasma.hpp"

int main( int argc, char** argv )
{
	std::map<std::string, double> parameterMap {
		{ "CentralCellField", 1.0 },
		{ "Throatfield", 10.0 },
		{ "PlasmaRadiusMin", 0.2 },
		{ "PlasmaRadiusMax", 0.6 }
	};
	auto pVacConfig = std::make_shared<MirrorPlasma::VacuumMirrorConfiguration>( parameterMap,FuelName,false,false,false,false,false,false,"","" );
	auto pReferencePlasmaState = std::make_shared<MirrorPlasma>(pVacuumConfig, parameterMap, VoltageTrace);

	std::cout << protonImpactIonizationCrossSection( 1.0 ) << std::endl;

	return 0;

}
