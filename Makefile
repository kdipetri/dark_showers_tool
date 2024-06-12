# Set the shell.
SHELL=/usr/bin/env bash
ROOT=`root-config --cflags --glibs`

# Include the configuration.
-include Makefile.inc

CXX_COMMON:=-I$(PREFIX_INCLUDE) $(CXX_COMMON) $(GZIP_LIB)
CXX_COMMON+= -L$(PREFIX_LIB) -Wl,-rpath,$(PREFIX_LIB) -lpythia8 -ldl
PYTHIA=$(PREFIX_LIB)/libpythia8$(LIB_SUFFIX)

# Rules without physical targets (secondary expansion for specific rules).
.SECONDEXPANSION:

# make all
all: mkdirBin card_runner

# make executable bin
mkdirBin:
	mkdir -p $(PWD)/bin

# verify PYTHIA library is installed
$(PYTHIA):
	$(error Error: PYTHIA must be built, please run in the top PYTHIA directory)

# higgs portal $(PYTHIA)
card_runner: src/$$@.cpp
	$(CXX) src/$@.cpp -o bin/$@.exe -w $(CXX_COMMON) $(ROOT_LIB) -ldl $(ROOT) $(HEPMC3_INCLUDE) $(HEPMC3_LIB)

# clean working dir
clean:
	rm -f *~
	rm -f \#*.*#
	rm -f $(PWD)/include/#*.*#
	rm -f $(PWD)/include/*~
	rm -f $(PWD)/src/#*.*#
	rm -f $(PWD)/src/*~
	rm -f $(PWD)/bin/*.exe
	rmdir bin
