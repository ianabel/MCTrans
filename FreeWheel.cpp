#include "SundialsWrapper.hpp" 


/*
 *
 * When doing a freewheeling plasma the time-dependent variables are 
 * V ( or equivalently the Mach number ), T_i, T_e
 *
 */

int ARKStep_FreeWheel( realtype t, N_Vector u, N_Vector uDot, void* voidPlasma )
{
	MirrorPlasma* plasmaPtr = reinterpret_cast<MirrorPlasma*>( voidPlasma );

	double TiOld = plasmaPtr->IonTemperature;
	double TeOld = plasmaPtr->ElectronTemperature;

	if ( ION_TEMPERATURE( u ) < 0.0 ) {
#if defined( DEBUG )
		std::cerr << "Error in SUNDIALS solve, due to negative ion temperature" << std::endl;
#endif
		return 1;
	}
	if ( ELECTRON_TEMPERATURE( u ) < 0.0 ) {
#if defined( DEBUG )
		std::cerr << "Error in SUNDIALS solve, due to negative electron temperature" << std::endl;
#endif
		return 2;
	}

	// We're now evolving the voltage with time 
	plasmaPtr->IonTemperature = ION_TEMPERATURE( u );
	plasmaPtr->ElectronTemperature = ELECTRON_TEMPERATURE( u );
	plasmaPtr->ImposedVoltage = VOLTAGE( u );


	
	try {
		plasmaPtr->SetTime(t);
	} catch ( std::domain_error &e ) {
		// Timestep too long?
#ifdef DEBUG
		std::cerr << "Evaluating RHS at t = " << std::setprecision( 20 ) << t << " ?!" << std::endl;
#endif
		return 3;
	}
	plasmaPtr->SetMachFromVoltage();
	plasmaPtr->ComputeSteadyStateNeutrals();
#if defined( DEBUG ) && defined( SUNDIALS_DEBUG ) && defined( INTERNAL_RK_DEBUG )
	std::cerr << "t = " << t << " ; T_i = " << plasmaPtr->IonTemperature << " ; T_e = " << plasmaPtr->ElectronTemperature << " MachNumber " << plasmaPtr->MachNumber << std::endl;
#endif

	// I d omega / dt = ( Change in Angular Momentum )
	// I d omega / dt = ( MomentumToVoltage )^(-1) dV/dt
	// so dV/dt = MtoV * ( Change in Angular Momentum )
	
	double MomentumToVoltage = plasmaPtr->PlasmaColumnWidth * plasmaPtr->PlasmaCentralRadius() * plasmaPtr->CentralCellFieldStrength / plasmaPtr->MomentOfInertia();

	try {
		double IonHeating  = plasmaPtr->IonHeating();
		double IonHeatLoss = plasmaPtr->IonHeatLosses();
		double ElectronHeating  = plasmaPtr->ElectronHeating();
		double ElectronHeatLoss = plasmaPtr->ElectronHeatLosses();

#if defined( DEBUG ) && defined( SUNDIALS_DEBUG ) && defined( INTERNAL_RK_DEBUG )
		std::cerr << " Ion Heating      = " << IonHeating      << " ; Ion Heat Loss       = " << IonHeatLoss      << std::endl;
		std::cerr << " Electron Heating = " << ElectronHeating << " ; Electron Heat Loss  = " << ElectronHeatLoss << std::endl;
#endif

		ION_HEAT_BALANCE( uDot )      = ( IonHeating - IonHeatLoss );
		ELECTRON_HEAT_BALANCE( uDot ) = ( ElectronHeating - ElectronHeatLoss );

		// Should be negative to decelerate the plasma
		double RadialCurrent = -VOLTAGE( u ) / plasmaPtr->ExternalResistance;
		double AngularMomentumInjection = plasmaPtr->InjectedTorque( RadialCurrent );
		double AngularMomentumLoss = plasmaPtr->TotalAngularMomentumLosses();

		MOMENTUM_BALANCE( uDot ) = ( MomentumToVoltage ) * ( AngularMomentumInjection - AngularMomentumLoss );


	} catch ( std::exception& e ) {
		return -1;
	} 

	plasmaPtr->IonTemperature = TiOld;
	plasmaPtr->ElectronTemperature = TeOld;
	plasmaPtr->SetMachFromVoltage();
	plasmaPtr->ComputeSteadyStateNeutrals();

	return ARK_SUCCESS;
}

// Let the plasma decay through a resistor
void MCTransConfig::doFreeWheel( MirrorPlasma& plasma ) const
{
	sundials::Context sunctx;	
	sunindextype NDims = 3; // Two temps & a voltage
	N_Vector initialCondition = N_VNew_Serial( NDims, sunctx );

	ION_TEMPERATURE( initialCondition ) = plasma.IonTemperature;
	ELECTRON_TEMPERATURE( initialCondition ) = plasma.ElectronTemperature;
	VOLTAGE( initialCondition ) = plasma.ImposedVoltage;

	realtype t0 = plasma.time;

	void *arkMem = ARKStepCreate( nullptr, ARKStep_FreeWheel, t0, initialCondition, sunctx );

	if ( arkMem == nullptr ) {
		throw std::runtime_error( "Cannot allocate ARKStep Working Memory" );
	}

	// Dummy Jacobian, will be filled by ARKStep with finite-difference approximations
	SUNMatrix       Jacobian = SUNDenseMatrix( NDims, NDims, sunctx );
	// Small system, direct solve is fastest
	SUNLinearSolver  LS = SUNLinSol_Dense( initialCondition, Jacobian, sunctx );

	ArkodeErrorWrapper( ARKStepSetLinearSolver( arkMem, LS, Jacobian ), "ARKStepSetLinearSolver" );
	
	

	double abstol = plasma.SundialsAbsTol;
	double reltol = plasma.SundialsRelTol;

#ifdef DEBUG
	std::cerr << "Using SundialsAbsTol = " << abstol << " and SundialsRelTol = " << reltol << std::endl;
#endif
	ArkodeErrorWrapper( ARKStepSStolerances( arkMem, reltol, abstol ), "ARKStepSStolerances" );
	ArkodeErrorWrapper( ARKStepSetTableNum( arkMem, IRK_SCHEME, ARKSTEP_NULL_STEPPER ), "ARKStepSetTableNum" );
	
	ArkodeErrorWrapper( ARKStepSetUserData( arkMem, reinterpret_cast<void*>( &plasma ) ), "ARKStepSetUserData" );


	const unsigned long MaxSteps = 1e4;
	ArkodeErrorWrapper( ARKStepSetMaxNumSteps( arkMem, MaxSteps ), "ARKStepSetMaxNumSteps" );

	realtype t,tRet = 0;	
	int errorFlag;

#ifdef DEBUG
	std::cerr << "Solving from t = 0 to t = " << EndTime << std::endl;
	std::cerr << "Writing output every " << OutputDeltaT << std::endl;
#endif 
	ArkodeErrorWrapper( ARKStepSetStopTime( arkMem, EndTime ), "ARKStepSetStopTime" );
	for ( t = t0 + OutputDeltaT; t < EndTime; t += OutputDeltaT )
	{
		errorFlag = ARKStepEvolve( arkMem, t, initialCondition, &tRet, ARK_NORMAL );
		switch ( errorFlag ) {
			case ARK_SUCCESS:
				break;
			default:
				throw std::runtime_error( "ARKStep failed with error " + std::to_string( errorFlag ) );
			break;
		}

		plasma.ElectronTemperature = ELECTRON_TEMPERATURE( initialCondition );
		plasma.IonTemperature = ION_TEMPERATURE( initialCondition );
		plasma.SetMachFromVoltage();
		plasma.ComputeSteadyStateNeutrals();
		plasma.WriteTimeslice( tRet );
		plasma.SetTime( tRet );
#if defined( DEBUG )
		std::cerr << "Writing timeslice at t = " << t << std::endl;
#endif
#if defined( DEBUG ) && defined( SUNDIALS_DEBUG )
	std::cerr << "After evolving to " << tRet << " T_i = " << ION_TEMPERATURE( initialCondition ) << " ; T_e = " << ELECTRON_TEMPERATURE( initialCondition ) << std::endl;
#endif
	}

	// We've solved and found the answer. Update the plasma object

	plasma.ElectronTemperature = ELECTRON_TEMPERATURE( initialCondition );
	plasma.IonTemperature      =      ION_TEMPERATURE( initialCondition );
	plasma.ImposedVoltage      =              VOLTAGE( initialCondition );

#ifdef DEBUG
	long nSteps = 0,nfeEvals = 0,nfiEvals = 0;
	ArkodeErrorWrapper( ARKStepGetNumSteps( arkMem, &nSteps ), "ARKGetNumSteps" );
	ArkodeErrorWrapper( ARKStepGetNumRhsEvals( arkMem, &nfeEvals, &nfiEvals ), "ARKGetNumRhsEvals" );
	std::cerr << "SUNDIALS Timestepping took " << nSteps << " internal timesteps resulting in " << nfiEvals << " implicit function evaluations" << std::endl;
#endif 
	// Teardown 
	{
		SUNLinSolFree( LS );
		SUNMatDestroy( Jacobian );
		N_VDestroy( initialCondition );
		ARKStepFree( &arkMem );
	}

	plasma.SetMachFromVoltage();
	plasma.ComputeSteadyStateNeutrals();
}
