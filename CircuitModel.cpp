#include "SundialsWrapper.hpp" 


/*
 *
 * When doing a freewheeling plasma the time-dependent variables are 
 * V ( or equivalently the Mach number ), T_i, T_e
 *
 */

#define N_DIMS_CIRCUITMODEL 5 /* I, V_cap, T_i, T_e, V_plasma */

#define V_CAP_IDX 3
#define I_CAP_IDX 4
#define V_CAP( u ) NV_Ith_S( u, V_CAP_IDX )
#define I_CAP( u ) NV_Ith_S( u, I_CAP_IDX )
#define V_CAP_EQN( F ) NV_Ith_S( F, V_CAP_IDX )
#define I_CAP_EQN( F ) NV_Ith_S( F, I_CAP_IDX )



int ARKStep_CircuitModel( realtype t, N_Vector u, N_Vector uDot, void* voidPlasma )
{
	MirrorPlasma* plasmaPtr = reinterpret_cast<MirrorPlasma*>( voidPlasma );

	// External circuit parameters.
	double C_cap = plasmaPtr->CBCapacitance;
	double R_cap = plasmaPtr->CBInternalResistance;
	double L_line = plasmaPtr->CBLineInductance;
	double R_line = plasmaPtr->CBLineResistance;

	double TiOld = plasmaPtr->IonTemperature;
	double TeOld = plasmaPtr->ElectronTemperature;
	double VOld = plasmaPtr->ImposedVoltage;

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
	plasmaPtr->UpdatePhi();
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
		double AngularMomentumInjection = plasmaPtr->InjectedTorque( I_CAP( u ) );
		double AngularMomentumLoss = plasmaPtr->TotalAngularMomentumLosses();

		MOMENTUM_BALANCE( uDot ) = ( MomentumToVoltage ) * ( AngularMomentumInjection - AngularMomentumLoss );

		V_CAP_EQN( uDot ) = -I_CAP( u )/C_cap - V_CAP( u )/( R_cap * C_cap );
		I_CAP_EQN( uDot ) = ( V_CAP( u ) - VOLTAGE( u ) - R_line * I_CAP( u ) ) / L_line;


	} catch ( std::exception& e ) {
		return -1;
	} 

	plasmaPtr->IonTemperature = TiOld;
	plasmaPtr->ElectronTemperature = TeOld;
	plasmaPtr->ImposedVoltage = VOld;
	plasmaPtr->SetMachFromVoltage();
	plasmaPtr->UpdatePhi();
	plasmaPtr->ComputeSteadyStateNeutrals();

	return ARK_SUCCESS;
}

// Let the plasma decay through a resistor
void MCTransConfig::doCircuitModel( MirrorPlasma& plasma ) const
{
	sundials::Context sunctx;	
	sunindextype NDims = N_DIMS_CIRCUITMODEL;
	N_Vector initialCondition = N_VNew_Serial( NDims, sunctx );

	ION_TEMPERATURE( initialCondition ) = plasma.IonTemperature;
	ELECTRON_TEMPERATURE( initialCondition ) = plasma.ElectronTemperature;
	VOLTAGE( initialCondition ) = plasma.ImposedVoltage;
	V_CAP( initialCondition ) = plasma.CBChargedVoltage;
	I_CAP( initialCondition ) = 0.0; 

	realtype t0 = plasma.time;

	void *arkMem = ARKStepCreate( nullptr, ARKStep_CircuitModel, t0, initialCondition, sunctx );

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

	N_Vector positivityEnforcement = N_VNew_Serial( NDims, sunctx );
	N_VConst( 0.0, positivityEnforcement ); // Default to no constraints
	ION_TEMPERATURE( positivityEnforcement ) = 2.0;      // T_i > 0
	ELECTRON_TEMPERATURE( positivityEnforcement ) = 2.0; // T_e > 0

	ArkodeErrorWrapper( ARKStepSetConstraints( arkMem, positivityEnforcement ), "ARKStepSetConstraints" );

	// Because the scheme is 4th order, we request cubic hermite interpolation between
	// internal timesteps, and don't allow the timestep to exceed 5*dt where dt is the
	// time between outputs.

	ArkodeErrorWrapper( ARKStepSetInterpolantDegree( arkMem, 3 ), "ARKStepSetInterpolantDegree" );
	ArkodeErrorWrapper( ARKStepSetMaxStep( arkMem, OutputDeltaT*5 ), "ARKStepSetMaxStep" );

	const unsigned long MaxSteps = 1e4;
	ArkodeErrorWrapper( ARKStepSetMaxNumSteps( arkMem, MaxSteps ), "ARKStepSetMaxNumSteps" );

	realtype t,tRet = 0;	
	int errorFlag;

#ifdef DEBUG
	std::cerr << "Solving from t = " << plasma.time << " to t = " << EndTime << std::endl;
	std::cerr << "Writing output every " << OutputDeltaT << std::endl;
#endif 
	ArkodeErrorWrapper( ARKStepSetStopTime( arkMem, EndTime ), "ARKStepSetStopTime" );
	for ( t = t0 + OutputDeltaT; t < EndTime; t += OutputDeltaT )
	{
#if defined( DEBUG )
		double curTime;
		ArkodeErrorWrapper( ARKStepGetCurrentTime( arkMem, &curTime ), "ARKStepGetCurrentTime" );
#endif
		if ( t > EndTime )
			t = EndTime;
		errorFlag = ARKStepEvolve( arkMem, t, initialCondition, &tRet, ARK_NORMAL );
		switch ( errorFlag ) {
			case ARK_SUCCESS:
#if defined( DEBUG )
				std::cerr << "Internal time is " << curTime << " Evolved to " << tRet << " with intent of reaching " << t << std::endl;
#endif
				break;
			default:
				throw std::runtime_error( "ARKStep failed with error " + std::to_string( errorFlag ) );
			break;
		}

		// ARKStep has evolved us to t = tRet, update the plasma object and write it out.
		plasma.SetTime( tRet );
		plasma.ElectronTemperature = ELECTRON_TEMPERATURE( initialCondition );
		plasma.IonTemperature = ION_TEMPERATURE( initialCondition );
		plasma.ImposedVoltage = VOLTAGE( initialCondition );
		plasma.SetMachFromVoltage();
		plasma.ComputeSteadyStateNeutrals();
		plasma.WriteTimeslice( tRet );
#if defined( DEBUG )
		std::cerr << "Writing timeslice at t = " << tRet << std::endl;
#endif
#if defined( DEBUG ) && defined( SUNDIALS_DEBUG )
		std::cerr << "After evolving to " << tRet << " T_i = " << ION_TEMPERATURE( initialCondition ) << " ; T_e = " << ELECTRON_TEMPERATURE( initialCondition ) << std::endl;
#endif
		double RelativeIonRate = ::fabs( ( plasma.IonHeating() - plasma.IonHeatLosses() )/( plasma.IonDensity * plasma.IonTemperature * ReferenceTemperature * ReferenceDensity ) );
		double RelativeElectronRate =::fabs( ( plasma.ElectronHeating() - plasma.ElectronHeatLosses() )/( plasma.ElectronDensity * plasma.ElectronTemperature * ReferenceTemperature * ReferenceDensity ) );
		if ( !plasma.isTimeDependent &&
		     RelativeIonRate < plasma.RateThreshold &&
		     RelativeElectronRate < plasma.RateThreshold )
		{
#if defined( DEBUG )
	std::cerr << "Steady state reached at time " << tRet << " with T_i = " << ION_TEMPERATURE( initialCondition ) << " ; T_e = " << ELECTRON_TEMPERATURE( initialCondition ) << std::endl;
#endif
			break;
		}
	}

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
}
