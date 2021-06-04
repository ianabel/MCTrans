
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

#define ION_TEMPERATURE( u ) NV_Ith_S( u, ION_TEMP_IDX )
#define ELECTRON_TEMPERATURE( u ) NV_Ith_S( u, ELEC_TEMP_IDX )
#define DENSITY( u ) NV_Ith_S( u, DENSITY_IDX )

#define ION_HEAT_BALANCE( F ) NV_Ith_S( F, ION_TEMP_IDX )
#define ELECTRON_HEAT_BALANCE( F ) NV_Ith_S( F, ELEC_TEMP_IDX )
#define PARTICLE_BALANCE( F ) NV_Ith_S( F, DENSITY_IDX )

#include "Config.hpp"
#include "MirrorPlasma.hpp"
#include <exception>
#include <iostream>

int KINSysWrapper( N_Vector u, N_Vector F, void* voidPlasma )
{
	MirrorPlasma* plasmaPtr = reinterpret_cast<MirrorPlasma*>( voidPlasma );

	plasmaPtr->IonTemperature = ION_TEMPERATURE( u );
	plasmaPtr->ElectronTemperature = ELECTRON_TEMPERATURE( u );
	//plasma->ElectronDensity = DENSITY( u );
	//plasma->SetIonDensity() // Set n_i from n_e, Z_i
	
	plasmaPtr->SetMachFromVoltage();

	std::cerr << "Iteration at T_i = " << plasmaPtr->IonTemperature << " ; T_e = " << plasmaPtr->ElectronTemperature << std::endl;


	try {
		double IonHeating  = plasmaPtr->IonHeating();
		double IonHeatLoss = plasmaPtr->IonHeatLosses();
		double ElectronHeating  = plasmaPtr->ElectronHeating();
		double ElectronHeatLoss = plasmaPtr->ElectronHeatLosses();

		std::cerr << " Ion Heating      = " << IonHeating << " ; Ion Heat Loss       = " << IonHeatLoss << std::endl;
		std::cerr << " Electron Heating = " << IonHeating << " ; Electron Heat Loss  = " << IonHeatLoss << std::endl;

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

	double InitialTemperature = 0.1; // Assume Ti = Te = 100eV to begin
	ION_TEMPERATURE( initialCondition ) = InitialTemperature;
	ELECTRON_TEMPERATURE( initialCondition ) = InitialTemperature;

	errorFlag = KINInit( kinMem, KINSysWrapper, initialCondition );

	// Dummy Jacobian, will be filled by KINSol with finite-difference approximations
	SUNMatrix       Jacobian = SUNDenseMatrix( NDims, NDims );
	// Small system, direct solve is fastest
	SUNLinearSolver  LS = SUNLinSol_Dense( initialCondition, Jacobian );

	errorFlag = KINSetLinearSolver( kinMem, LS, Jacobian );

	if ( errorFlag != KIN_SUCCESS ) {
		throw std::runtime_error( "Could not initialise KINSol" );
	}

	double ftol = 1.e-6;
	double scstol = 1.e-6;
	double jtol = 1.e-6;

	errorFlag = KINSetFuncNormTol( kinMem, ftol );
	errorFlag = KINSetScaledStepTol( kinMem, scstol );
	errorFlag = KINSetRelErrFunc( kinMem, jtol );

	N_Vector positivityEnforcement = N_VNew_Serial( NDims );
	N_VConst( 0.0, positivityEnforcement ); // Default to no constraints
	ION_TEMPERATURE( positivityEnforcement ) = 1.0;      // T_i >= 0
	ELECTRON_TEMPERATURE( positivityEnforcement ) = 1.0; // T_e >= 0

	errorFlag = KINSetConstraints( kinMem, positivityEnforcement );

	errorFlag = KINSetUserData( kinMem, reinterpret_cast<void*>( &plasma ) );

	// set one to be the constant vector containing N 1s
	// this turns off scaling
	N_Vector one = N_VNew_Serial( NDims );
	N_VConst( 1.0, one );

	errorFlag = KINSetMaxSetupCalls( kinMem, 1 );
	errorFlag = KINSetPrintLevel( kinMem, 1 );
	errorFlag = KINSol( kinMem, initialCondition, KIN_LINESEARCH, one, one );

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

