FC=gfortran
FFLAGS=-O2 -std=f2018 -Wall -fimplicit-none -Jbuild

# NOTE: fpm is the primary build system and handles module dependency order
# automatically.  This Makefile is a best-effort fallback only.
# The wildcard SRC list below does NOT guarantee a dependency-correct compile
# order.  As inter-module dependencies grow, replace the wildcard with an
# explicit, dependency-ordered list of source files.
SRC=$(wildcard src/*.f90)

build:
	mkdir -p build
	$(FC) $(FFLAGS) -c $(SRC)
	mv *.o build/ 2>/dev/null || true

clean:
	rm -rf build *.mod *.o
