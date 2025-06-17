FC = gfortran
FLAGS = -O2 -fopenmp
BESSEL_UTILS_DIR= ../tests
OBJS = bessel_functions_toms.o bessel_utils.o io_utils.o WH_solver.o 
TARGET = WH_solver

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(FLAGS) -o $@ $(OBJS)

bessel_functions_toms.o: $(BESSEL_UTILS_DIR)/bessel_functions_toms.f90
	$(FC) $(FLAGS) -c $(BESSEL_UTILS_DIR)/bessel_functions_toms.f90

bessel_utils.o:  $(BESSEL_UTILS_DIR)/bessel_utils.f90 bessel_functions_toms.o
	$(FC) $(FLAGS) -c  $(BESSEL_UTILS_DIR)/bessel_utils.f90

io_utils.o:  $(BESSEL_UTILS_DIR)/io_utils.f90 
	$(FC) $(FLAGS) -c  $(BESSEL_UTILS_DIR)/io_utils.f90

WH_solver.o: WH_solver.f90 bessel_utils.o io_utils.o
	$(FC) $(FLAGS) -c WH_solver.f90

clean:
	rm -f *.o *.mod $(TARGET)

