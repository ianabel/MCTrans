
all: MCTrans++

include Makefile.config


SOURCES = MCTransConfig.cpp MirrorPlasma.cpp FusionYield.cpp Report.cpp AlphaHeating.cpp Neutrals.cpp SundialsWrapper.cpp NetCDFIO.cpp BatchRunner.cpp FreeWheel.cpp
HEADERS = MirrorPlasma.hpp FusionYield.hpp Config.hpp Species.hpp PlasmaPhysics.hpp NetCDFIO.hpp AtomicPhysics.hpp BatchRunner.hpp SundialsWrapper.hpp
OBJECTS = $(patsubst %.cpp,%.o,$(SOURCES))


%.o: %.cpp Makefile Makefile.local $(HEADERS)
	$(CXX) -c $(CXXFLAGS) -o $@ $<

MCTrans++: MCTrans.o $(OBJECTS) $(HEADERS) Makefile Makefile.local
	$(CXX) $(CXXFLAGS) -o MCTrans++ MCTrans.o $(OBJECTS) $(SUN_LINK_FLAGS) $(NETCDF_LINK_FLAGS)

MCTrans.pdf: manual/Makefile manual/MCTrans.tex manual/macros.tex manual/references.bib
	make -C manual MCTrans.pdf
	ln -s manual/MCTrans.pdf

test: MCTrans++
	cd examples/; ./check_examples

unit_tests: unit_test_suite
	./unit_test_suite

TEST_SOURCES = tests/MCTransTests.cpp tests/NeutralTests.cpp

CXX_TEST_FLAGS = $(CXXFLAGS)

unit_test_suite: $(TEST_SOURCES)	$(OBJECTS) $(HEADERS) Makefile Makefile.local
	$(CXX) $(CXX_TEST_FLAGS) -o unit_test_suite $(TEST_SOURCES) $(OBJECTS) $(SUN_LINK_FLAGS) $(NETCDF_LINK_FLAGS)

clean: 
	rm -f MCTrans++ $(OBJECTS)

distclean: clean
	rm -f MCTrans.pdf

.PHONY: examples clean all test distclean
.SUFFIXES:
