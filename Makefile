
all: MCTrans++

MCTrans++: MCTrans.cpp Centrifugal.cpp Centrifugal.hpp Config.hpp FusionYield.hpp FusionYield.cpp AlphaHeating.cpp AlphaHeating.hpp
	g++ -g -Og --std=c++17 -Wall -o MCTrans++ MCTrans.cpp Centrifugal.cpp FusionYield.cpp AlphaHeating.cpp
