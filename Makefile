
CXX ?= g++
CXXFLAGS += -g -Og --std=c++20 -Wall
all: MCTrans++ MCTrans.pdf

MCTrans++: MCTrans.cpp MirrorPlasma.hpp MirrorPlasma.cpp FusionYield.hpp FusionYield.cpp Report.cpp AlphaHeating.cpp Config.hpp PlasmaPhysics.hpp  Species.hpp
	$(CXX) $(CXXFLAGS) -o MCTrans++ MCTrans.cpp MirrorPlasma.cpp FusionYield.cpp AlphaHeating.cpp Report.cpp 

MCTrans.pdf: manual/MCTrans.tex manual/macros.tex manual/references.bib
	make -C manual MCTrans.pdf
	ln -s manual/MCTrans.pdf

examples: examples/*.report

examples/%.report: examples/%.conf MCTrans++
	./MCTrans++ $< > $@

clean: 
	rm MCTrans++ MCTrans.pdf

.PHONY: examples clean all
.SUFFIXES:
