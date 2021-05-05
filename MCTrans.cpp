#include <iostream>
#include <iomanip>
#include <cmath>
#include <memory>
#include <boost/math/tools/roots.hpp>

#include "Config.hpp"
#include "MirrorPlasma.hpp"


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

	MCTransConfig config( fname );

	std::unique_ptr<MirrorPlasma> result = config.Solve();
	result->PrintReport();

	return 0;
}

std::unique_ptr<MirrorPlasma> MCTransConfig::Solve() const
{
	std::unique_ptr<MirrorPlasma> plasma = std::make_unique<MirrorPlasma>( *ReferencePlasmaState );
	switch ( Type ) {
		case SolveType::SteadyStateMachSolve:
			doMachSolve( *plasma );
			break;
		case SolveType::SteadyStateTempSolve:
			doTempSolve( *plasma );
			break;
		default:
			throw std::invalid_argument( "Unknown Solve Type" );
	}
	
	return plasma;
}

void MCTransConfig::doMachSolve( MirrorPlasma& plasma ) const
{
	// NB This uses power densities in W/m^3
	auto PowerBalance = [ &plasma ]( double M ) {
		plasma.MachNumber = M;

		double HeatLoss = plasma.IonHeatLosses() + plasma.ElectronHeatLosses();
		double Heating = plasma.IonHeating() + plasma.ElectronHeating();

		return Heating - HeatLoss;
	};

	boost::uintmax_t iters = 1000;
	boost::math::tools::eps_tolerance<double> tol( 11 ); // only bother getting part in 1024 accuracy
	double InitialMach = 1.0; // Usually M > 4 for these solutions
	double Factor = 2.0;
	bool rising = true; // Confinement gets uniformly better for increasing M, and Viscous heating increases with M
	auto [ M_lower, M_upper ] = boost::math::tools::bracket_and_solve_root( PowerBalance, InitialMach, Factor, rising, tol, iters );
	plasma.MachNumber = ( M_lower + M_upper )/2.0;
}

void MCTransConfig::doTempSolve( MirrorPlasma& plasma ) const
{
	const double TiTe = plasma.IonTemperature / plasma.ElectronTemperature;
	// NB This uses power densities in W/m^3
	auto PowerBalance = [ &plasma, TiTe ]( double Te ) {
		plasma.ElectronTemperature = Te;
		plasma.IonTemperature = Te * TiTe;

		// Update Mach Number from new T_e
		plasma.SetMachFromVoltage();
		
		double HeatLoss = plasma.IonHeatLosses() + plasma.ElectronHeatLosses();
		double Heating = plasma.IonHeating() + plasma.ElectronHeating();

		return Heating - HeatLoss;
	};

	boost::uintmax_t iters = 1000;
	boost::math::tools::eps_tolerance<double> tol( 11 ); // only bother getting part in 1024 accuracy
	double InitialTe = 0.1; // Start cold -- 100eV
	double Factor = 2.0;
	bool rising = false; // At fixed Mach Number, heating the plasma up will increase losses
	auto [ T_lower, T_upper ] = boost::math::tools::bracket_and_solve_root( PowerBalance, InitialTe, Factor, rising, tol, iters );
	plasma.ElectronTemperature = ( T_lower + T_upper )/2.0;
	plasma.SetMachFromVoltage();
}
