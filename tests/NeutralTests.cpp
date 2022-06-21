
// No BOOST_TEST_MODULE, that's in MCTransTests
#include <boost/test/unit_test.hpp>

#include "../MirrorPlasma.hpp"
#include "../AtomicPhysics.hpp"


BOOST_AUTO_TEST_SUITE( atomic_physics_test_suite, * boost::unit_test::tolerance( 1e-6 ) )


BOOST_AUTO_TEST_CASE( electron_impact_ionization_test )
{
	// Check it's 0 below the min energy
	BOOST_TEST( electronImpactIonizationCrossSection( 13.5 ) == 0.0 );
	// Check zero energy gives zero rather than an error
	BOOST_TEST( electronImpactIonizationCrossSection( 0.0 ) == 0.0 );
	// Check -ve energy gives zero rather than an error
	BOOST_TEST( electronImpactIonizationCrossSection( -1.0 ) == 0.0 );

	std::vector<std::pair<double,double>> JanevReferenceData{
		{ 1.50E+01, 7.70E-18 },
		{ 2.00E+01, 2.96E-17 },
		{ 4.00E+01, 5.88E-17 },
		{ 6.00E+01, 6.19E-17 },
		{ 8.00E+01, 5.88E-17 },
		{ 1.00E+02, 5.51E-17 },
		{ 2.00E+02, 3.95E-17 },
		{ 4.00E+02, 2.41E-17 },
		{ 6.00E+02, 1.75E-17 },
		{ 8.00E+02, 1.39E-17 },
		{ 1.00E+03, 1.15E-17 },
		{ 2.00E+03, 6.26E-18 },
		{ 4.00E+03, 3.51E-18 },
		{ 6.00E+03, 2.43E-18 },
		{ 8.00E+03, 1.88E-18 },
		{ 1.00E+04, 1.53E-18 }
	};

	for ( auto sample : JanevReferenceData )
		BOOST_TEST( electronImpactIonizationCrossSection( sample.first ) == sample.second, boost::test_tools::tolerance( 0.035 ) );
}



BOOST_AUTO_TEST_SUITE_END()
