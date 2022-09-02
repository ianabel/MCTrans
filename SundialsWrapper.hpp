
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

void ArkodeErrorWrapper( int errorFlag, std::string&& fName );

