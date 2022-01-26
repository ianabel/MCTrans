
all: MCTrans++

include Makefile.config

SOURCES = MCTrans.cpp MirrorPlasma.cpp FusionYield.cpp Report.cpp AlphaHeating.cpp Neutrals.cpp SundialsWrapper.cpp NetCDFIO.cpp
HEADERS = MirrorPlasma.hpp FusionYield.hpp Config.hpp Species.hpp PlasmaPhysics.hpp NetCDFIO.hpp
OBJECTS = $(patsubst %.cpp,%.o,$(SOURCES))


%.o: %.cpp Makefile Makefile.local $(HEADERS)
	$(CXX) -c $(CXXFLAGS) -o $@ $<

MCTrans++: $(OBJECTS) $(HEADERS) Makefile Makefile.local
	$(CXX) $(CXXFLAGS) -o MCTrans++ $(OBJECTS) $(SUN_LINK_FLAGS) $(NETCDF_LINK_FLAGS)

MCTrans.pdf: manual/Makefile manual/MCTrans.tex manual/macros.tex manual/references.bib
	make -C manual MCTrans.pdf
	ln -s manual/MCTrans.pdf

examples: examples/*.report


examples/%.report: examples/%.conf MCTrans++
	$(warning "WARNING: Replacing the old output $@ with current output of MCTrans++ $<")
	./MCTrans++ $< > $@

test: MCTrans++
	cd examples/; ./check_examples

clean: 
	rm -f MCTrans++ $(OBJECTS)

distclean: clean
	rm -f MCTrans.pdf

.PHONY: examples clean all test distclean
.SUFFIXES:
