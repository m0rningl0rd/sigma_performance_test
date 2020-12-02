#!make

FFLAGS = -O3 -mkl=sequential -fp-model fast=2 -no-prec-div -static-intel
# FFLAGS = -O2 -mkl=sequential -fp-model precise -static-intel

F90 = ifort

sigma_performance_test: sigma_performance_test.o
	${F90} ${FFLAGS} sigma_performance_test.o -o sigma_performance_test

sigma_performance_test.o: sigma_performance_test.f90
	${F90} -c ${FFLAGS} sigma_performance_test.f90

clean:
	@rm -f *.o
