
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
}

BOOST_AUTO_TEST_SUITE_END()
