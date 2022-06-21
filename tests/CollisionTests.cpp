

#include <algorithm>
#include <complex>
#include <vector>
#include <utility>

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include <boost/math/constants/constants.hpp>

// This defines the name of the test suite and causes
// the default behaviour of the BOOST_TEST macros to compare
// within a part in 10^6 rather than exact comparison.
BOOST_AUTO_TEST_SUITE( collisional_transport_test_suite, * boost::unit_test::tolerance( 1e-6 ) )

// Compare our function for the coulomb logarithm against a set of known values
// computed from the formulary.
BOOST_AUTO_TEST_CASE( log_lambda_test )
{
	
}

BOOST_AUTO_TEST_SUITE_END()
