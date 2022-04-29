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

std::shared_ptr<MirrorPlasma> MCTransConfig::Solve()
{
	switch ( Type ) {
		case SolveType::SteadyStateMachSolve:
			doMachSolve( *ReferencePlasmaState );
			break;
		case SolveType::SteadyStateTempSolve:
			doTempSolve( *ReferencePlasmaState );
			break;
		default:
			throw std::invalid_argument( "Unknown Solve Type" );
	}

	return std::move( ReferencePlasmaState );
}

void MCTransConfig::doMachSolve( MirrorPlasma& plasma ) const
{
	// NB This uses power densities in W/m^3
	auto PowerBalance = [ &plasma ]( double M ) {
		plasma.MachNumber = M;
		plasma.ComputeSteadyStateNeutrals();

		double HeatLoss = plasma.IonHeatLosses() + plasma.ElectronHeatLosses();
		double Heating = plasma.IonHeating() + plasma.ElectronHeating();

		return Heating - HeatLoss;
	};

	boost::uintmax_t iters = 1000;
	boost::math::tools::eps_tolerance<double> tol( 11 ); // only bother getting part in 1024 accuracy
	double InitialMach = plasma.initialMach(); // Usually M > 4 for these solutions
	double Factor = 1.25;
	bool rising = true; // Confinement gets uniformly better for increasing M, and Viscous heating increases with M
	auto M_bounds = boost::math::tools::bracket_and_solve_root( PowerBalance, InitialMach, Factor, rising, tol, iters );
	auto M_lower = M_bounds.first, M_upper = M_bounds.second;
	plasma.MachNumber = ( M_lower + M_upper )/2.0;
	plasma.ComputeSteadyStateNeutrals();
}


