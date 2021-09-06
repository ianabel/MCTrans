
include Makefile.local

CXX ?= g++

STD=c++17

ifdef DEBUG
CXXFLAGS += -DDEBUG -g3 -O0 --std=$(STD) -Werror -pedantic
else
CXXFLAGS += -O2 --std=$(STD) -Wall
endif

all: MCTrans++

SUNDIALS_DIR ?= /usr/local

SUNDIALS_INC=$(SUNDIALS_DIR)/include
SUNDIALS_LIB=$(SUNDIALS_DIR)/lib

SUNFLAGS=-I$(SUNDIALS_INC)
SUN_LINK_FLAGS = -L$(SUNDIALS_LIB) -Wl,-rpath $(SUNDIALS_LIB) -lsundials_arkode -lsundials_nvecserial

TOML11_DIR ?= ./toml11
TOML_FLAGS = -I$(TOML11_DIR)

CXXFLAGS += $(TOML_FLAGS) $(SUNFLAGS)

ifdef BOOST_DIR
	BOOST_FLAGS = -I$(BOOST_DIR)
	CXXFLAGS += $(BOOST_FLAGS)
endif

NETCDF_LINK_FLAGS =

ifdef NETCDF_DIR
	CXXFLAGS += -I$(NETCDF_DIR)/include
	NETCDF_LINK_FLAGS = -L$(NETCDF_DIR)/lib -Wl,-rpath $(NETCDF_DIR)/lib
endif

ifdef NETCDF_CXX_DIR
	CXXFLAGS += -I$(NETCDF_CXX_DIR)/include
	NETCDF_LINK_FLAGS += -L$(NETCDF_CXX_DIR)/lib -Wl,-rpath $(NETCDF_CXX_DIR)/lib
endif
	

NETCDF_LINK_FLAGS += -lnetcdf -lnetcdf_c++4

SOURCES = MCTrans.cpp MirrorPlasma.cpp FusionYield.cpp Report.cpp AlphaHeating.cpp Neutrals.cpp SundialsWrapper.cpp NetCDFIO.cpp
HEADERS = MirrorPlasma.hpp FusionYield.hpp Config.hpp Species.hpp PlasmaPhysics.hpp NetCDFIO.hpp
OBJECTS = $(patsubst %.cpp,%.o,$(SOURCES))

Makefile.local:
	$(error You need to provide a Makefile.local for your machine. Try copying Makefile.local.example)

%.o: %.cpp Makefile Makefile.local $(HEADERS)
	$(CXX) -c $(CXXFLAGS) -o $@ $<

MCTrans++: $(OBJECTS) $(HEADERS) Makefile Makefile.local
	$(CXX) $(CXXFLAGS) -o MCTrans++ $(OBJECTS) $(SUN_LINK_FLAGS) $(NETCDF_LINK_FLAGS)

MCTrans.pdf: manual/Makefile manual/MCTrans.tex manual/macros.tex manual/references.bib
	make -C manual MCTrans.pdf
	ln -s manual/MCTrans.pdf

examples: examples/*.report


examples/%.report: examples/%.conf MCTrans++
	$(warning "WARNING: Replacing the old output $@ with current output of MCTrans++ $<)
	./MCTrans++ $< > $@

test: MCTrans++
	cd examples/; ./check_examples

clean: 
	rm -f MCTrans++ $(OBJECTS)

distclean: clean
	rm -f MCTrans.pdf

.PHONY: examples clean all test distclean
.SUFFIXES:
