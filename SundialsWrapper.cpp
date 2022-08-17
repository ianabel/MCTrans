
#include <arkode/arkode_arkstep.h>     /* access to ARKode func., consts. */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector       */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix       */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype */

// Number of equations/unknowns in system
#define N_DIMENSIONS 2 

#define ION_TEMP_IDX 0
#define ELEC_TEMP_IDX 1
#define DENSITY_IDX 2
#define MOMENTUM_IDX 2 

#define MACH_IDX 1

#define ION_TEMPERATURE( u ) NV_Ith_S( u, ION_TEMP_IDX )
#define ELECTRON_TEMPERATURE( u ) NV_Ith_S( u, ELEC_TEMP_IDX )
#define DENSITY( u ) NV_Ith_S( u, DENSITY_IDX )
#define VOLTAGE( u ) NV_Ith_S( u, MOMENTUM_IDX )

#define MACH_NUMBER( u ) NV_Ith_S( u, MACH_IDX )

#define ION_HEAT_BALANCE( F ) NV_Ith_S( F, ION_TEMP_IDX )
#define ELECTRON_HEAT_BALANCE( F ) NV_Ith_S( F, ELEC_TEMP_IDX )
#define PARTICLE_BALANCE( F ) NV_Ith_S( F, DENSITY_IDX )
#define MOMENTUM_BALANCE( F ) NV_Ith_S( F, MOMENTUM_IDX )

#define IRK_SCHEME static_cast<ARKODE_DIRKTableID>( ARKSTEP_DEFAULT_DIRK_4 )
#define ARKSTEP_NULL_STEPPER  static_cast<ARKODE_ERKTableID>( -1 )

#include "Config.hpp"
#include "MirrorPlasma.hpp"
#include <exception>
#include <iostream>
#include <iomanip>

void ArkodeErrorWrapper( int errorFlag, std::string&& fName )
{
	if ( errorFlag == ARK_SUCCESS )
		return;
	else
	{
		std::string errName = ARKStepGetReturnFlagName( errorFlag );
		throw std::runtime_error( "Error " + errName + " returned from ARKode function: " + fName );
	}
}

int ARKStep_TemperatureSolve( realtype t, N_Vector u, N_Vector uDot, void* voidPlasma )
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

	plasmaPtr->IonTemperature = ION_TEMPERATURE( u );
	plasmaPtr->ElectronTemperature = ELECTRON_TEMPERATURE( u );
	//plasma->ElectronDensity = DENSITY( u );
	//plasma->SetIonDensity() // Set n_i from n_e, Z_i
	
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
//		PARTICLE_BALANCE( uDot ) = ParticleBalance; 


	} catch ( std::exception& e ) {
		return -1;
	} 

	plasmaPtr->IonTemperature = TiOld;
	plasmaPtr->ElectronTemperature = TeOld;
	plasmaPtr->SetMachFromVoltage();
	plasmaPtr->UpdatePhi();
	plasmaPtr->ComputeSteadyStateNeutrals();

	return ARK_SUCCESS;
}

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

// In this mode, Mach Number is u(0) and T_i is u(1)
// but we still solve both the power-balance equations by running to steady state
int ARKStep_FixedTeSolve( realtype t, N_Vector u, N_Vector F, void* voidPlasma )
{
	MirrorPlasma* plasmaPtr = reinterpret_cast<MirrorPlasma*>( voidPlasma );

	plasmaPtr->MachNumber = MACH_NUMBER( u );

	plasmaPtr->IonTemperature = ION_TEMPERATURE( u );
	

	try {
		double IonHeating  = plasmaPtr->IonHeating();
		double IonHeatLoss = plasmaPtr->IonHeatLosses();
		double ElectronHeating  = plasmaPtr->ElectronHeating();
		double ElectronHeatLoss = plasmaPtr->ElectronHeatLosses();

		ION_HEAT_BALANCE( F )      = ( IonHeating - IonHeatLoss );
		ELECTRON_HEAT_BALANCE( F ) = ( ElectronHeating - ElectronHeatLoss );
//		PARTICLE_BALANCE( F ) = ParticleBalance; 

	} catch ( std::exception& e ) {
		return -1;
	} 

	return 0;
}


void MCTransConfig::doTempSolve( MirrorPlasma& plasma ) const
{
	sundials::Context sunctx;
	sunindextype NDims = N_DIMENSIONS;
	N_Vector initialCondition = N_VNew_Serial( NDims, sunctx );

	// Set initial temperature to be such that the mach number is reasonable.
	double InitialMach = plasma.InitialMach;
	double SoundSpeed = ::fabs( plasma.ImposedVoltage / ( plasma.PlasmaColumnWidth * plasma.CentralCellFieldStrength * InitialMach ) );
	double InitialTemperature = plasma.IonSpecies.Mass * ProtonMass * SoundSpeed * SoundSpeed / ReferenceTemperature;
#ifdef DEBUG
	std::cerr << "Setting T_0 = " << InitialTemperature  << " to have M at t=0 fixed to " << InitialMach << std::endl;
#endif


	ION_TEMPERATURE( initialCondition ) = InitialTemperature;
	ELECTRON_TEMPERATURE( initialCondition ) = InitialTemperature;
	
	plasma.ElectronTemperature = InitialTemperature;
	plasma.IonTemperature = InitialTemperature;
	plasma.SetMachFromVoltage();
	plasma.UpdatePhi();
	plasma.ComputeSteadyStateNeutrals();

	plasma.InitialiseNetCDF();

	realtype t0 = 0;

	void *arkMem = ARKStepCreate( nullptr, ARKStep_TemperatureSolve, t0, initialCondition, sunctx );

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

	ArkodeErrorWrapper( ARKStepSetMaxStep( arkMem, OutputDeltaT*5 ), "ARKStepSetMaxStep" );

	const unsigned long MaxSteps = 1e5;
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
#if defined( DEBUG )
		double curTime;
		ArkodeErrorWrapper( ARKStepGetCurrentTime( arkMem, &curTime ), "ARKStepGetCurrentTime" );
#endif
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

		plasma.SetTime( tRet );
		plasma.ElectronTemperature = ELECTRON_TEMPERATURE( initialCondition );
		plasma.IonTemperature = ION_TEMPERATURE( initialCondition );
		plasma.SetMachFromVoltage();
		plasma.UpdatePhi();
		plasma.ComputeSteadyStateNeutrals();
		plasma.WriteTimeslice( tRet );
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
	
	plasma.FinaliseNetCDF();

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

// e.g. Feedback Control of Voltage on Te
void MCTransConfig::doFixedTeSolve( MirrorPlasma& plasma ) const
{
	throw std::logic_error( "Fixed T_e solve not yet fully implemented" );
	sundials::Context sunctx;
	sunindextype NDims = N_DIMENSIONS;
	N_Vector initialCondition = N_VNew_Serial( NDims, sunctx );

	double InitialTemperature = plasma.InitialTemp;
	ION_TEMPERATURE( initialCondition ) = InitialTemperature;
	ELECTRON_TEMPERATURE( initialCondition ) = InitialTemperature;

	realtype t0 = 0;

	void *arkMem = ARKStepCreate( nullptr, ARKStep_TemperatureSolve, t0, initialCondition, sunctx );

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

	ArkodeErrorWrapper( ARKStepSStolerances( arkMem, reltol, abstol ), "ARKStepSStolerances" );
	ArkodeErrorWrapper( ARKStepSetTableNum( arkMem, IRK_SCHEME, ARKSTEP_NULL_STEPPER ), "ARKStepSetTableNum" );
	
	ArkodeErrorWrapper( ARKStepSetUserData( arkMem, reinterpret_cast<void*>( &plasma ) ), "ARKStepSetUserData" );


	const unsigned long MaxSteps = 1e4;
	ArkodeErrorWrapper( ARKStepSetMaxNumSteps( arkMem, MaxSteps ), "ARKStepSetMaxNumSteps" );

	realtype t,tRet;	
	int errorFlag = ARKStepEvolve( arkMem, t, initialCondition, &tRet, ARK_NORMAL );
	switch ( errorFlag ) {
		case ARK_SUCCESS:
			break;
		default:
			throw std::runtime_error( "KINSol failed with error " + std::to_string( errorFlag ) );
			break;
	}	

	// We've solved and found the answer. Update the plasma object

	plasma.ElectronTemperature = ELECTRON_TEMPERATURE( initialCondition );
	plasma.IonTemperature      =      ION_TEMPERATURE( initialCondition );

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

// Let the plasma decay through a resistor
void MCTransConfig::doFreeWheel( MirrorPlasma& plasma ) const
{
	sundials::Context sunctx;	
	sunindextype NDims = N_DIMENSIONS;
	N_Vector initialCondition = N_VNew_Serial( NDims, sunctx );

	double InitialTemperature = plasma.InitialTemp;
	ION_TEMPERATURE( initialCondition ) = InitialTemperature;
	ELECTRON_TEMPERATURE( initialCondition ) = InitialTemperature;

	plasma.ElectronTemperature = InitialTemperature;
	plasma.IonTemperature = InitialTemperature;
	plasma.MachNumber = plasma.InitialMach;

	plasma.ImposedVoltage = plasma.MachNumber * plasma.SoundSpeed() * plasma.PlasmaColumnWidth * plasma.CentralCellFieldStrength;

	plasma.ComputeSteadyStateNeutrals();

	plasma.InitialiseNetCDF();

	realtype t0 = 0;

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
	}

	// We've solved and found the answer. Update the plasma object

	plasma.ElectronTemperature = ELECTRON_TEMPERATURE( initialCondition );
	plasma.IonTemperature      =      ION_TEMPERATURE( initialCondition );
	plasma.SetTime( t );

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
