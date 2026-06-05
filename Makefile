FC=gfortran
FFLAGS=-O2 -std=f2008 -Wall -fimplicit-none -Jbuild
SRC=$(wildcard src/*.f90)
build:
	mkdir -p build
	$(FC) $(FFLAGS) -c $(SRC)
	mv *.o build/ 2>/dev/null || true
clean:
	rm -rf build *.mod *.o
