
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

SUNFLAGS=-I$(SUNDIALS_INC)
SUN_LINK_FLAGS = -L$(SUNDIALS_LIB) -Wl,-rpath $(SUNDIALS_LIB) -lsundials_kinsol -lsundials_nvecserial

TOML11_DIR ?= ./toml11
TOML_FLAGS = -I$(TOML11_DIR)

CXXFLAGS += $(TOML_FLAGS) $(SUNFLAGS)

SOURCES = MCTrans.cpp MirrorPlasma.cpp FusionYield.cpp Report.cpp AlphaHeating.cpp Neutrals.cpp SundialsWrapper.cpp
HEADERS = MirrorPlasma.hpp FusionYield.hpp Config.hpp Species.hpp PlasmaPhysics.hpp
OBJECTS = $(patsubst %.cpp,%.o,$(SOURCES))

Makefile.local:
	$(error You need to provide a Makefile.local for your machine. Try copying Makefile.local.example)

%.o: %.cpp Makefile Makefile.local $(HEADERS)
	$(CXX) -c $(CXXFLAGS) -o $@ $<

MCTrans++: $(OBJECTS) $(HEADERS) Makefile Makefile.local
	$(CXX) $(CXXFLAGS) -o MCTrans++ $(OBJECTS) $(SUN_LINK_FLAGS)

MCTrans.pdf: manual/Makefile manual/MCTrans.tex manual/macros.tex manual/references.bib
	make -C manual MCTrans.pdf
	ln -s manual/MCTrans.pdf

examples: examples/*.report

examples/%.report: examples/%.conf MCTrans++
	./MCTrans++ $< > $@

test: MCTrans++
	true

clean: 
	rm -f MCTrans++ $(OBJECTS)

distclean: clean
	rm -f MCTrans.pdf

.PHONY: examples clean all test distclean
.SUFFIXES:
