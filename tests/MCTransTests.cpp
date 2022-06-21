
#define BOOST_TEST_MODULE MCTrans
#include <boost/test/included/unit_test.hpp>

#include "../MirrorPlasma.hpp"

// Defines the entry point and a simple set of functional tests for the MirrorPlasma class

// This defines the name of the test suite and causes
// the default behaviour of the BOOST_TEST macros to compare
// within a part in 10^6 rather than exact comparison.
BOOST_AUTO_TEST_SUITE( functiontal_test_suite, * boost::unit_test::tolerance( 1e-6 ) );

BOOST_AUTO_TEST_CASE( construct_mirror_plasma_test )
{
	std::map<std::string, double> parameters{
		{ "CentralCellField", 1.0 }
	};
	
	MirrorPlasma *pSamplePlasma = nullptr;

	BOOST_REQUIRE_NO_THROW( pSamplePlasma = new MirrorPlasma( parameters, "Hydrogen", false, false, false, false, false, false, "", "", "" ) );

	BOOST_TEST( pSamplePlasma->CentralCellFieldStrength == 1.0 );

}

BOOST_AUTO_TEST_SUITE_END();
