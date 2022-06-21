
void MCTransConfig::FreeWheelSolve( MirrorPlasma& plasma ) const
{
	
	sunindextype NDims = N_DIMENSIONS;
	N_Vector initialCondition = N_VNew_Serial( NDims );

	double InitialTemperature = plasma.InitialTemp;
	ION_TEMPERATURE( initialCondition ) = InitialTemperature;
	ELECTRON_TEMPERATURE( initialCondition ) = InitialTemperature;

	plasma.ElectronTemperature = InitialTemperature;
	plasma.IonTemperature = InitialTemperature;
	plasma.SetMachFromVoltage();
	plasma.ComputeSteadyStateNeutrals();

	plasma.InitialiseNetCDF();

	realtype t0 = 0;

	void *arkMem = ARKStepCreate( nullptr, ARKStep_TemperatureSolve, t0, initialCondition );

	if ( arkMem == nullptr ) {
		throw std::runtime_error( "Cannot allocate ARKStep Working Memory" );
	}

	// Dummy Jacobian, will be filled by ARKStep with finite-difference approximations
	SUNMatrix       Jacobian = SUNDenseMatrix( NDims, NDims );
	// Small system, direct solve is fastest
	SUNLinearSolver  LS = SUNLinSol_Dense( initialCondition, Jacobian );

	ArkodeErrorWrapper( ARKStepSetLinearSolver( arkMem, LS, Jacobian ), "ARKStepSetLinearSolver" );
	
	

	double abstol = plasma.SundialsAbsTol;
	double reltol = plasma.SundialsRelTol;

#ifdef DEBUG
	std::cerr << "Using SundialsAbsTol = " << abstol << " and SundialsRelTol = " << reltol << std::endl;
#endif
	ArkodeErrorWrapper( ARKStepSStolerances( arkMem, reltol, abstol ), "ARKStepSStolerances" );
	ArkodeErrorWrapper( ARKStepSetTableNum( arkMem, DEFAULT_DIRK_5, -1 ), "ARKStepSetTableNum" );
	
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
	for ( t = OutputDeltaT; t < EndTime; t += OutputDeltaT )
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
		plasma.WriteTimeslice( t );
#if defined( DEBUG )
		std::cerr << "Writing timeslice at t = " << t << std::endl;
#endif
#if defined( DEBUG ) && defined( SUNDIALS_DEBUG )
	std::cerr << "After evolving to " << tRet << " T_i = " << ION_TEMPERATURE( initialCondition ) << " ; T_e = " << ELECTRON_TEMPERATURE( initialCondition ) << std::endl;
#endif

		double RelativeIonRate = ::fabs( ( plasma.IonHeating() - plasma.IonHeatLosses() )/( plasma.IonDensity * plasma.IonTemperature * ReferenceTemperature * ReferenceDensity ) );
		double RelativeElectronRate =::fabs( ( plasma.ElectronHeating() - plasma.ElectronHeatLosses() )/( plasma.ElectronDensity * plasma.ElectronTemperature * ReferenceTemperature * ReferenceDensity ) );
#if defined( DEBUG ) && defined( SUNDIALS_DEBUG )
		std::cerr << " Relative Rate of Change in Ion Energy Density " << RelativeIonRate * 100 << " %/s" << std::endl;
		std::cerr << " Relative Rate of Change in Electron Energy Density " << RelativeElectronRate * 100 << " %/s" << std::endl;
#endif
		if ( !plasma.isTimeDependent &&
		     RelativeIonRate < plasma.RateThreshold &&
		     RelativeElectronRate < plasma.RateThreshold )
		{
#if defined( DEBUG )
	std::cerr << "Steady state reached at time " << t << " with T_i = " << ION_TEMPERATURE( initialCondition ) << " ; T_e = " << ELECTRON_TEMPERATURE( initialCondition ) << std::endl;
#endif
			break;
		}

	}

	// We've solved and found the answer. Update the plasma object

	plasma.ElectronTemperature = ELECTRON_TEMPERATURE( initialCondition );
	plasma.IonTemperature      =      ION_TEMPERATURE( initialCondition );
	// plasma.SetTime( EndTime );

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
