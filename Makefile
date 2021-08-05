
include Makefile.local

CXX ?= g++

STD=c++17

ifdef DEBUG
CXXFLAGS += -DDEBUG -g -Og --std=$(STD) -Werror -pedantic
else
CXXFLAGS += -O2 --std=$(STD) -Wall
endif

all: MCTrans++

SUNDIALS_DIR ?= /usr/local

SUNDIALS_INC=$(SUNDIALS_DIR)/include
SUNDIALS_LIB=$(SUNDIALS_DIR)/lib

SUNFLAGS=-I$(SUNDIALS_INC) -L$(SUNDIALS_LIB) -Wl,-rpath $(SUNDIALS_LIB) 
SUN_LINK_FLAGS = -lsundials_kinsol -lsundials_nvecserial

TOML11_DIR ?= ./toml11
TOML_FLAGS = -I$(TOML11_DIR)

CXXFLAGS += $(TOML_FLAGS) $(SUNFLAGS)

OBJECTS = MCTrans.cpp MirrorPlasma.cpp FusionYield.cpp Report.cpp AlphaHeating.cpp Neutrals.cpp SundialsWrapper.cpp
HEADERS = MirrorPlasma.hpp FusionYield.hpp Config.hpp Species.hpp PlasmaPhysics.hpp

MCTrans++: $(OBJECTS) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o MCTrans++ $(OBJECTS) $(SUN_LINK_FLAGS)

MCTrans.pdf: manual/MCTrans.tex manual/macros.tex manual/references.bib
	make -C manual MCTrans.pdf
	ln -s manual/MCTrans.pdf

examples: examples/*.report

examples/%.report: examples/%.conf MCTrans++
	./MCTrans++ $< > $@

test: MCTrans++
	true

clean: 
	rm -f MCTrans++

distclean: clean
	rm -f MCTrans.pdf

.PHONY: examples clean all test
.SUFFIXES:
