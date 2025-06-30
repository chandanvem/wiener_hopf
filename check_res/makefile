FC = gfortran
FLAGS = -O3 -fopenmp
BESSEL_UTILS_DIR= ../tests
OBJS = bessel_functions_toms.o bessel_utils.o io_utils.o check_residue.o
TARGET = check_residue

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(FLAGS) -o $@ $(OBJS)

bessel_functions_toms.o: $(BESSEL_UTILS_DIR)/bessel_functions_toms.f90
	$(FC) $(FLAGS) -c $(BESSEL_UTILS_DIR)/bessel_functions_toms.f90

bessel_utils.o:  $(BESSEL_UTILS_DIR)/bessel_utils.f90 bessel_functions_toms.o
	$(FC) $(FLAGS) -c  $(BESSEL_UTILS_DIR)/bessel_utils.f90

io_utils.o:  $(BESSEL_UTILS_DIR)/io_utils.f90
	$(FC) $(FLAGS) -c  $(BESSEL_UTILS_DIR)/io_utils.f90

check_residue.o: check_residue.f90 bessel_utils.o io_utils.o
	$(FC) $(FLAGS) -c check_residue.f90

clean:
	rm -f *.o *.mod $(TARGET)

