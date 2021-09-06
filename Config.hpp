#ifndef CONFIG_HPP
#define CONFIG_HPP

#include "MirrorPlasma.hpp"
#include <toml.hpp>
#include <memory>

class MCTransConfig {
	public:
		MCTransConfig( std::string const& configFile ) {
			const toml::value config = toml::parse( configFile );

			ReferencePlasmaState = std::make_unique<MirrorPlasma>( config );

			bool MachSolve = ( ReferencePlasmaState->ElectronTemperature > 0.0 );
			bool TempSolve = ( ReferencePlasmaState->pVacuumConfig->ImposedVoltage > 0.0 );
			if ( MachSolve && !TempSolve )
				Type = SteadyStateMachSolve;
			else if ( !MachSolve && TempSolve )
				Type = SteadyStateTempSolve;
			else if ( MachSolve && TempSolve )
				throw std::invalid_argument( "[error] Cannot specify both temperature and voltage." );
			else if ( !MachSolve && !TempSolve )
				throw std::invalid_argument( "[error] Must specify at least one of ElectronTemperature or Voltage." );

			if ( config.count( "timestepping" ) ) {
				const auto & timestep_conf = toml::find<toml::table>( config, "timestepping" );
				OutputDeltaT = timestep_conf.at( "OutputCadence" ).as_floating();
				EndTime = timestep_conf.at( "EndTime" ).as_floating();
			} else {
				OutputDeltaT = 0.0005; // 500µs cadence
				EndTime = 0.100; // 100 ms run time
			}
		};
		enum SolveType {
			SteadyStateMachSolve,
			SteadyStateTempSolve
		} Type;

		std::unique_ptr<MirrorPlasma> Solve();
		std::unique_ptr<MirrorPlasma> ReferencePlasmaState;
	private:
		void doMachSolve( MirrorPlasma& plasma ) const;
		void doTempSolve( MirrorPlasma& plasma ) const;
		void doFixedTeSolve( MirrorPlasma& plasma ) const;

		double OutputDeltaT,EndTime;

};

#endif // CONFIG_HPP
