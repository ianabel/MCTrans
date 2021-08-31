
#include <kinsol/kinsol.h>             /* access to KINSOL func., consts. */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector       */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix       */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype */

// Number of equations/unknowns in system
#define N_DIMENSIONS 2 

#define ION_TEMP_IDX 0
#define ELEC_TEMP_IDX 1
#define DENSITY_IDX 2

#define MACH_IDX 0 

#define ION_TEMPERATURE( u ) NV_Ith_S( u, ION_TEMP_IDX )
#define ELECTRON_TEMPERATURE( u ) NV_Ith_S( u, ELEC_TEMP_IDX )
#define DENSITY( u ) NV_Ith_S( u, DENSITY_IDX )

#define MACH_NUMBER( u ) NV_Ith_S( u, MACH_IDX )

#define ION_HEAT_BALANCE( F ) NV_Ith_S( F, ION_TEMP_IDX )
#define ELECTRON_HEAT_BALANCE( F ) NV_Ith_S( F, ELEC_TEMP_IDX )
#define PARTICLE_BALANCE( F ) NV_Ith_S( F, DENSITY_IDX )

#include "Config.hpp"
#include "MirrorPlasma.hpp"
#include <exception>
#include <iostream>

void KinsolErrorWrapper( int errorFlag, std::string&& fName )
{
	if ( errorFlag == KIN_SUCCESS )
		return;
	else
	{
		throw std::runtime_error( "Error " + std::to_string( errorFlag ) + " returned from KINSol function: " + fName );
	}
}



int KINSysWrapper_TemperatureSolve( N_Vector u, N_Vector F, void* voidPlasma )
{
	MirrorPlasma* plasmaPtr = reinterpret_cast<MirrorPlasma*>( voidPlasma );

	plasmaPtr->IonTemperature = ION_TEMPERATURE( u );
	plasmaPtr->ElectronTemperature = ELECTRON_TEMPERATURE( u );
	//plasma->ElectronDensity = DENSITY( u );
	//plasma->SetIonDensity() // Set n_i from n_e, Z_i
	
	plasmaPtr->SetMachFromVoltage();
#if defined( DEBUG ) && defined( SUNDIALS_DEBUG )
	std::cerr << "Iteration at T_i = " << plasmaPtr->IonTemperature << " ; T_e = " << plasmaPtr->ElectronTemperature << " MachNumber " << plasmaPtr->MachNumber << std::endl;
#endif


	try {
		double IonHeating  = plasmaPtr->IonHeating();
		double IonHeatLoss = plasmaPtr->IonHeatLosses();
		double ElectronHeating  = plasmaPtr->ElectronHeating();
		double ElectronHeatLoss = plasmaPtr->ElectronHeatLosses();

#if defined( DEBUG ) && defined( SUNDIALS_DEBUG )
		std::cerr << " Ion Heating      = " << IonHeating      << " ; Ion Heat Loss       = " << IonHeatLoss      << std::endl;
		std::cerr << " Electron Heating = " << ElectronHeating << " ; Electron Heat Loss  = " << ElectronHeatLoss << std::endl;
#endif

		ION_HEAT_BALANCE( F )      = ( IonHeating - IonHeatLoss );
		ELECTRON_HEAT_BALANCE( F ) = ( ElectronHeating - ElectronHeatLoss );
//		PARTICLE_BALANCE( F ) = ParticleBalance; 

	} catch ( std::exception& e ) {
		return -1;
	} 

	return 0;
}

// In this mode, Mach Number is u(0) and T_i is u(1)
// but we still solve both the power-balance equations
int KINSysWrapper_FixedTeSolve( N_Vector u, N_Vector F, void* voidPlasma )
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
	void *kinMem = KINCreate();
	if ( kinMem == nullptr ) {
		throw std::runtime_error( "Cannot allocate KINSol Working Memory" );
	}

	sunindextype NDims = N_DIMENSIONS;
	int errorFlag = KIN_SUCCESS;

	N_Vector initialCondition = N_VNew_Serial( NDims );

	double InitialTemperature = plasma.pVacuumConfig->InitialTemp;
	ION_TEMPERATURE( initialCondition ) = InitialTemperature;
	ELECTRON_TEMPERATURE( initialCondition ) = InitialTemperature;

	errorFlag = KINInit( kinMem, KINSysWrapper_TemperatureSolve, initialCondition );

	if ( errorFlag != KIN_SUCCESS ) {
		throw std::runtime_error( "Could not initialise KINSol" );
	}

	// Dummy Jacobian, will be filled by KINSol with finite-difference approximations
	SUNMatrix       Jacobian = SUNDenseMatrix( NDims, NDims );
	// Small system, direct solve is fastest
	SUNLinearSolver  LS = SUNLinSol_Dense( initialCondition, Jacobian );

	KinsolErrorWrapper( KINSetLinearSolver( kinMem, LS, Jacobian ), "KINSetLinearSolver" );
	
	

	double ftol = 1.e-3;
	double scstol = 1.e-11;
	double jtol = 1.e-9;

	KinsolErrorWrapper( KINSetFuncNormTol( kinMem, ftol ), "KINSetFuncNormTol" );
	KinsolErrorWrapper( KINSetScaledStepTol( kinMem, scstol ), "KINSetScaledStepTol" );
	KinsolErrorWrapper( KINSetRelErrFunc( kinMem, jtol ), "KINSetRelErrFunc" );

	N_Vector positivityEnforcement = N_VNew_Serial( NDims );
	N_VConst( 0.0, positivityEnforcement ); // Default to no constraints
	ION_TEMPERATURE( positivityEnforcement ) = 1.0;      // T_i >= 0
	ELECTRON_TEMPERATURE( positivityEnforcement ) = 1.0; // T_e >= 0

	KinsolErrorWrapper( KINSetConstraints( kinMem, positivityEnforcement ), "KINSetConstraints" );

	KinsolErrorWrapper( KINSetUserData( kinMem, reinterpret_cast<void*>( &plasma ) ), "KINSetUserData( kinMem, reinterpret_cast<void*>" );

	// set one to be the constant vector containing N 1s
	// this turns off scaling
	N_Vector one = N_VNew_Serial( NDims );
	N_VConst( 1.0, one );

	KinsolErrorWrapper( KINSetMaxSetupCalls( kinMem, 1 ), "KINSetMaxSetupCalls" );
#if defined( DEBUG ) && defined( SUNDIALS_DEBUG )
	KinsolErrorWrapper( KINSetPrintLevel( kinMem, 1 ), "KINSetPrintLevel" );
#else
	KinsolErrorWrapper( KINSetPrintLevel( kinMem, 0 ), "KINSetPrintLevel" );
#endif
	
	// Reference value for F(u) should be roughly
	// 1.5 * Final Temperature * Final Density / tau, where tau is some typical equilibrium timescale
	// choosing tau = 100ms, T = 1keV, n = 10^20 we have
	N_Vector f_scale = N_VNew_Serial( NDims );
	ELECTRON_TEMPERATURE( f_scale ) = 10 / ( 1.5 * ReferenceTemperature * ReferenceDensity );
	ION_TEMPERATURE( f_scale ) = 10 / ( 1.5 * ReferenceTemperature * ReferenceDensity );

#if defined( DEBUG ) && defined( SUNDIALS_DEBUG )
	std::cerr << "f_scale set to (" << ION_TEMPERATURE( f_scale ) << ", " << ELECTRON_TEMPERATURE( f_scale ) << ")" << std::endl;
#endif


	errorFlag =  KINSol( kinMem, initialCondition, KIN_LINESEARCH, one, f_scale );
	switch ( errorFlag ) {
		case KIN_SUCCESS:
			break;
		case KIN_INITIAL_GUESS_OK:
#ifdef DEBUG
			std::cerr << " Initial Temperature is a steady-state solution?! " << std::endl;
#endif
			break;
		default:
			throw std::runtime_error( "KINSol failed with error " + std::to_string( errorFlag ) );
			break;
	}	// We've solved and found the answer. Update the plasma object
	plasma.ElectronTemperature = ELECTRON_TEMPERATURE( initialCondition );
	plasma.IonTemperature      =      ION_TEMPERATURE( initialCondition );

	// Teardown 
	{
		SUNLinSolFree( LS );
		SUNMatDestroy( Jacobian );
		N_VDestroy( one );
		N_VDestroy( initialCondition );
		KINFree( &kinMem );
	}

	plasma.SetMachFromVoltage();
	plasma.ComputeSteadyStateNeutrals();
}

// e.g. Feedback Control of Voltage on Te
void MCTransConfig::doFixedTeSolve( MirrorPlasma& plasma ) const
{
	void *kinMem = KINCreate();
	if ( kinMem == nullptr ) {
		throw std::runtime_error( "Cannot allocate KINSol Working Memory" );
	}

	sunindextype NDims = N_DIMENSIONS;
	int errorFlag = KIN_SUCCESS;

	N_Vector initialCondition = N_VNew_Serial( NDims );

	double InitialTemperature = 0.1;
	ION_TEMPERATURE( initialCondition ) = InitialTemperature;
	ELECTRON_TEMPERATURE( initialCondition ) = InitialTemperature;

	errorFlag = KINInit( kinMem, KINSysWrapper_TemperatureSolve, initialCondition );
	if ( errorFlag != KIN_SUCCESS ) {
		throw std::runtime_error( "Could not initialise KINSol" );
	}

	// Dummy Jacobian, will be filled by KINSol with finite-difference approximations
	SUNMatrix       Jacobian = SUNDenseMatrix( NDims, NDims );
	// Small system, direct solve is fastest
	SUNLinearSolver  LS = SUNLinSol_Dense( initialCondition, Jacobian );

	KinsolErrorWrapper( KINSetLinearSolver( kinMem, LS, Jacobian ), "KINSetLinearSolver" );

	

	double ftol = 1.e-3;
	double scstol = 1.e-9;
	double jtol = 1.e-9;

	KinsolErrorWrapper( KINSetFuncNormTol( kinMem, ftol ), "KINSetFuncNormTol" );
	KinsolErrorWrapper( KINSetScaledStepTol( kinMem, scstol ), "KINSetScaledStepTol" );
	KinsolErrorWrapper( KINSetRelErrFunc( kinMem, jtol ), "KINSetRelErrFunc" );

	N_Vector positivityEnforcement = N_VNew_Serial( NDims );
	N_VConst( 0.0, positivityEnforcement ); // Default to no constraints
	ION_TEMPERATURE( positivityEnforcement ) = 1.0;      // T_i >= 0
	ELECTRON_TEMPERATURE( positivityEnforcement ) = 1.0; // T_e >= 0

	KinsolErrorWrapper( KINSetConstraints( kinMem, positivityEnforcement ), "KINSetConstraints" );

	KinsolErrorWrapper( KINSetUserData( kinMem, reinterpret_cast<void*>( &plasma ) ), "KINSetUserData( kinMem, reinterpret_cast<void*>" );

	// set one to be the constant vector containing N 1s
	// this turns off scaling
	N_Vector one = N_VNew_Serial( NDims );
	N_VConst( 1.0, one );

	KinsolErrorWrapper( KINSetMaxSetupCalls( kinMem, 1 ), "KINSetMaxSetupCalls" );
	KinsolErrorWrapper( KINSetPrintLevel( kinMem, 0 ), "KINSetPrintLevel" );

	// Reference value for F(u) should be roughly
	// 1.5 * Final Temperature * Final Density / tau, where tau is some typical equilibrium timescale
	// choosing tau = 100ms, T = 1keV, n = 10^20 we have
	N_Vector f_scale = N_VNew_Serial( NDims );
	ELECTRON_TEMPERATURE( f_scale ) = 10 / ( 1.5 * ReferenceTemperature * ReferenceDensity );
	ION_TEMPERATURE( f_scale ) = 10 / ( 1.5 * ReferenceTemperature * ReferenceDensity );


	errorFlag =  KINSol( kinMem, initialCondition, KIN_LINESEARCH, one, f_scale );
	switch ( errorFlag ) {
		case KIN_SUCCESS:
			break;
		case KIN_INITIAL_GUESS_OK:
#ifdef DEBUG
			std::cerr << " Initial Temperature is a steady-state solution?! " << std::endl;
#endif
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
		N_VDestroy( one );
		N_VDestroy( initialCondition );
		KINFree( &kinMem );
	}

	plasma.SetMachFromVoltage();
	plasma.ComputeSteadyStateNeutrals();
}
