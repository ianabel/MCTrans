
CXX ?= g++
CXXFLAGS += -g -Og --std=c++17 -Wall
all: MCTrans++ MCTrans.pdf

MCTrans++: MCTrans.cpp Centrifugal.cpp Centrifugal.hpp Config.hpp FusionYield.hpp FusionYield.cpp AlphaHeating.cpp AlphaHeating.hpp
	$(CXX) $(CXXFLAGS) -o MCTrans++ MCTrans.cpp Centrifugal.cpp FusionYield.cpp AlphaHeating.cpp

MCTrans.pdf: manual/MCTrans.tex manual/macros.tex manual/references.bib
	make -C manual MCTrans.pdf
	ln -s manual/MCTrans.pdf

clean: 
	rm MCTrans++ MCTrans.pdf

.PHONY: clean all
