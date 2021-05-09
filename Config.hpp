#ifndef CONFIG_HPP
#define CONFIG_HPP

#include "MirrorPlasma.hpp"
#include "toml11/toml.hpp"
#include <memory>

class MCTransConfig {
	public:
		MCTransConfig( std::string const& configFile ) {
			const toml::value config = toml::parse( configFile );

			if ( config.count( "plasma" ) != 1 )
				throw std::invalid_argument( "There must be one, and only one, [plasma] block in the configuration file." );

			const auto plasmaConfig = toml::find( config, "plasma" );

			ReferencePlasmaState = std::make_unique<MirrorPlasma>( plasmaConfig );

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
		};
		enum SolveType {
			SteadyStateMachSolve,
			SteadyStateTempSolve
		} Type;

		std::unique_ptr<MirrorPlasma> Solve() const;
		std::unique_ptr<MirrorPlasma> ReferencePlasmaState;
	private:
		void doMachSolve( MirrorPlasma& plasma ) const;
		void doTempSolve( MirrorPlasma& plasma ) const;

};

#endif // CONFIG_HPP
