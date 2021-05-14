
CXX ?= g++
CXXFLAGS += -g -Og --std=c++20 -Wall
all: MCTrans++ MCTrans.pdf

OBJECTS = MCTrans.cpp MirrorPlasma.cpp FusionYield.cpp Report.cpp AlphaHeating.cpp Neutrals.cpp
HEADERS = MirrorPlasma.hpp FusionYield.hpp Config.hpp Species.hpp PlasmaPhysics.hpp

MCTrans++: $(OBJECTS) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o MCTrans++ $(OBJECTS)

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
