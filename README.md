# flare
Fast LIGO/LISA Analysis of Response and Estimation

Data used to set up the Reduced Order Model for EOBNRv2HM should be untared and put in a directory pointed to by the environment variable ROM_DATA_PATH.

Data representing the noise (square root) PSD for the LIGO/VIRGO detectors must be located in a directory pointed to by the environment variable LLV_NOISE_DATA_PATH.

## Compiling the code

Flare depends on [BAMBI](https://ccpforge.cse.rl.ac.uk/gf/project/skynet/)
which needs to be installed first. The BAMBI configure script might create
some headaches due to missing `cblas_dgemm`. On Debian this dependency is
provided by `libatlas3-base`. Otherwise one can install ATLAS manually to
a location `$X` and then configure BAMBI with something like this:
```bash
LD_LIBRARY_PATH=${X}/lib:$LD_LIBRARY_PATH \
CFLAGS="-I${X}/include -L${X}/lib" \
CPPFLAGS="-I${X}/include -L${X}/lib" \
CXXFLAGS="-I${X}/include -L${X}/lib" \
./configure --prefix=${X}
```

Next, create a new `MACHINE` tag at the beginning of the flare Makefile, and
a corresponding new section in the `ifeq`-`else`-`endif` block.  In the new
section, set the following variables appropriately:
```bash
MESSAGE="Compiling for my local configuration"
GSLROOT = /usr
BAMBIROOT = /usr/local
FFTWROOT = /usr
FC = gfortran
CC = gcc
CXX = g++
CPP = g++
CFLAGS += -g -O3 -fopenmp -I $(GSLROOT)/include -I $(FFTWROOT)/include -L$(FFTWROOT)/lib -L$(GSLROOT)/lib
CXXFLAGS = -g -O3 -fopenmp
LD = $(CPP)
LDFLAGS = -g -C -fopenmp -L$(FFTWROOT)/lib -L$(GSLROOT)/lib
PTMCMC=$(PWD)/ptmcmc
```
The provided example settings are for a standard Debian system with GSL and
FFTW installed via the package manager and BAMBI installed in `/usr/local`.

Run `make` in the base directory of the flare repository.
