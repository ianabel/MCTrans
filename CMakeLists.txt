cmake_minimum_required(VERSION 3.10)

project(MCTrans++ VERSION 0.5)

set(MCTRANS_SRCS MCTrans.cpp MirrorPlasma.cpp FusionYield.cpp Report.cpp AlphaHeating.cpp Neutrals.cpp SundialsWrapper.cpp NetCDFIO.cpp)
set(MCTRANS_HEADERS MirrorPlasma.hpp FusionYield.hpp Config.hpp Species.hpp PlasmaPhysics.hpp NetCDFIO.hpp)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(MCTrans++ ${MCTRANS_SRCS} ${MCTRANS_HEADERS})

target_link_libraries(MCTrans++ PUBLIC sundials_arkode sundials_nvecserial netcdf netcdf_c++4)



