
#include <memory>
#include <boost/math/tools/roots.hpp>

#include "Config.hpp"
#include "MirrorPlasma.hpp"

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


