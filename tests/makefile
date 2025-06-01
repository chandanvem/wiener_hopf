# Makefile for Bessel project

FC = gfortran
FLAGS = -O2

OBJS = bessel_functions_toms.o bessel_utils.o bessel_test.o
TARGET = bessel_program

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(FLAGS) -o $@ $(OBJS)

bessel_functions_toms.o: bessel_functions_toms.f90
	$(FC) $(FLAGS) -c bessel_functions_toms.f90

bessel_utils.o: bessel_utils.f90 bessel_functions_toms.o
	$(FC) $(FLAGS) -c bessel_utils.f90

bessel_test.o: bessel_test.f90 bessel_utils.o
	$(FC) $(FLAGS) -c bessel_test.f90

clean:
	rm -f *.o *.mod $(TARGET)

