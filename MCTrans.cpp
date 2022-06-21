#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <memory>
#include <boost/math/tools/roots.hpp>

#include "Config.hpp"
#include "MirrorPlasma.hpp"
#include "BatchRunner.hpp"

extern "C" {
#include <fenv.h>
}

int main( int argc, char** argv )
{
	std::string fname( "Mirror.conf" );
	if ( argc == 2 )
		fname = argv[ 1 ];
	if ( argc > 2 )
	{
		std::cerr << "Usage: MCTrans++ ConfigFile.conf" << std::endl;
		return 1;
	}

#if defined(DEBUG) && defined(FPEXCEPT)
	std::cerr << "Floating Point Exceptions Enabled" << std::endl;
	::feenableexcept( FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW );
#endif


	BatchRunner runner(fname);
	runner.runBatchSolve();
	return 0;

}

