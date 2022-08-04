
#include <memory>
#include <boost/math/tools/roots.hpp>

#include "Config.hpp"
#include "MirrorPlasma.hpp"

std::shared_ptr<MirrorPlasma> MCTransConfig::Solve()
{
	MirrorPlasma& plasma = *ReferencePlasmaState;



	switch ( Type ) {
		case SolveType::SteadyStateMachSolve:
			doMachSolve( plasma );
			break;
		case SolveType::SteadyStateTempSolve:
			InitialisePlasma();
			plasma.InitialiseNetCDF();
			doTempSolve( plasma );
			break;
		case SolveType::FreewheelSolve:
			// Set initial conditions
			InitialisePlasma();
			plasma.InitialiseNetCDF();
			if ( ReferencePlasmaState->ImposedVoltage > 0.0 ) {
				std::cerr << "Evolving to steady state before decaying" << std::endl;
				plasma.isTimeDependent = false; // Just get a steady state
				doTempSolve( plasma );
			} else {
				throw std::invalid_argument( "Currently the way Free-wheel decay works is to run to a steady state first, then decay. Please set ImposedVoltage for the quasi-steady phase of the simulation" );
			}
			// Let the plasma spin down
			plasma.isTimeDependent = true; // Now run for the prescribed time
			EndTime += plasma.time; // So the time runs from 0 -> steady-state-time -> s-s-t + Endtime, giving a full 'EndTime' of decay
			doFreeWheel( plasma );
			break;
		default:
			throw std::invalid_argument( "Unknown Solve Type" );
	}

	plasma.FinaliseNetCDF();
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

void MCTransConfig::InitialisePlasma() const
{
	MirrorPlasma& plasma = *ReferencePlasmaState;
	double InitialMach = plasma.initialMach();
	double InitialTemperature = plasma.initialTemperature();

	if ( plasma.ImposedVoltage > 0.0 ) {
		// Set initial temperature to be such that the mach number is reasonable.
		if ( InitialMach > 0.0 ) {
			double SoundSpeed = ::fabs( plasma.ImposedVoltage / ( plasma.PlasmaColumnWidth * plasma.CentralCellFieldStrength * InitialMach ) );
			InitialTemperature = plasma.IonSpecies.Mass * ProtonMass * SoundSpeed * SoundSpeed / ReferenceTemperature;
#ifdef DEBUG
			std::cerr << "Setting T_0 = " << InitialTemperature  << " to have M at t=0 fixed to " << InitialMach << std::endl;
#endif
		}
	} else {
		if ( InitialMach > 0.0 && InitialTemperature > 0.0 ) {
			plasma.ElectronTemperature = InitialTemperature;
			plasma.ImposedVoltage = plasma.PlasmaColumnWidth * InitialMach * plasma.SoundSpeed() * plasma.CentralCellFieldStrength;
		}
	}
	plasma.ElectronTemperature = InitialTemperature;
	plasma.IonTemperature = InitialTemperature;
	plasma.SetMachFromVoltage();
	plasma.ComputeSteadyStateNeutrals();
}


