
// No BOOST_TEST_MODULE, that's in MCTransTests
#include <boost/test/unit_test.hpp>

#include "../MirrorPlasma.hpp"
#include "../AtomicPhysics.hpp"


BOOST_AUTO_TEST_SUITE( atomic_physics_test_suite, * boost::unit_test::tolerance( 1e-6 ) )


/* 
 * Tests for all the cross-sections include testing against the Janev data, with varying tolerance
 * the Janev data is often not good to more than 10-20% so the fit should be required to be much more accurate
 */

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

BOOST_AUTO_TEST_CASE( proton_impact_ionization_test )
{
	// Check it's 0 below the min energy
	BOOST_TEST( protonImpactIonizationCrossSection( 450.0 ) == 0.0 );
	// Check zero energy gives zero rather than an error
	BOOST_TEST( protonImpactIonizationCrossSection( 0.0 ) == 0.0 );
	// Check -ve energy gives zero rather than an error
	BOOST_TEST( protonImpactIonizationCrossSection( -1.0 ) == 0.0 );

	std::vector<std::pair<double,double>> JanevReferenceData{
		{ 5.00E+02, 1.46E-20 },
		{ 1.00E+03, 1.46E-19 },
		{ 2.00E+03, 1.02E-18 },
		{ 5.00E+03, 6.24E-18 },
		{ 1.00E+04, 1.94E-17 },
		{ 2.00E+04, 6.73E-17 },
		{ 5.00E+04, 1.43E-16 },
		{ 1.00E+05, 1.10E-16 },
		{ 2.00E+05, 6.99E-17 },
		{ 5.00E+05, 3.48E-17 },
		{ 1.00E+06, 1.94E-17 },
		{ 2.00E+06, 1.05E-17 },
		{ 5.00E+06, 4.62E-18 }
	};

	for ( auto sample : JanevReferenceData )
		BOOST_TEST( protonImpactIonizationCrossSection( sample.first ) == sample.second, boost::test_tools::tolerance( 0.05 ) );
}

BOOST_AUTO_TEST_CASE( hydrogen_charge_exchange_test )
{
	// Check it's 0 below the min energy
	BOOST_TEST( HydrogenChargeExchangeCrossSection( 0.1 ) == 0.0 );
	// Check zero energy gives zero rather than an error
	BOOST_TEST( HydrogenChargeExchangeCrossSection( 0.0 ) == 0.0 );
	// Check -ve energy gives zero rather than an error
	BOOST_TEST( HydrogenChargeExchangeCrossSection( -1.0 ) == 0.0 );

	// We're only using n=1 cross section, so this is from 2.3.1 of Janev
	std::vector<std::pair<double,double>> JanevReferenceData{
		{ 1.20E-01, 4.96E-15 },
		{ 2.00E-01, 4.70E-15 },
		{ 5.00E-01, 4.33E-15 },
		{ 1.00E+00, 4.10E-15 },
		{ 2.00E+00, 3.83E-15 },
		{ 5.00E+00, 3.46E-15 },
		{ 1.00E+01, 3.17E-15 },
		{ 2.00E+01, 2.93E-15 },
		{ 5.00E+01, 2.65E-15 },
		{ 1.00E+02, 2.44E-15 },
		{ 2.00E+02, 2.22E-15 },
		{ 5.00E+02, 1.97E-15 },
		{ 1.00E+03, 1.71E-15 },
		{ 2.00E+03, 1.44E-15 },
		{ 5.00E+03, 1.10E-15 },
		{ 1.00E+04, 7.75E-16 },
		{ 2.00E+04, 4.45E-16 },
		{ 5.00E+04, 9.93E-17 },
		{ 1.00E+05, 1.01E-17 },
		{ 2.00E+05, 6.09E-19 },
		{ 5.00E+05, 6.03E-21 },
		{ 1.00E+06, 1.57E-22 },
		{ 2.00E+06, 3.78E-24 },
		{ 5.00E+06, 2.56E-26 },
		{ 1.00E+07, 5.99E-28 }
	};

	for ( auto sample : JanevReferenceData )
		BOOST_TEST( HydrogenChargeExchangeCrossSection( sample.first ) == sample.second, boost::test_tools::tolerance( 0.15 ) );
}
BOOST_AUTO_TEST_SUITE_END()
