FC = gfortran
FLAGS = -O3 -fopenmp
BESSEL_UTILS_DIR= ./bessel_utils
IO_UTILS_DIR=./io_utils
INPUT_PARAMS_DIR=./input_param_utils
CONTOUR_UTILS_DIR=./contour_utils
OBJS = bessel_functions_toms.o 	bessel_utils.o\
	 input_param_utils.o io_utils.o\
	contour_init_utils.o contour_generate_utils.o\
	main.o 

TARGET = main

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(FLAGS) -o $@ $(OBJS)

bessel_functions_toms.o: $(BESSEL_UTILS_DIR)/bessel_functions_toms.f90
	$(FC) $(FLAGS) -c $(BESSEL_UTILS_DIR)/bessel_functions_toms.f90

bessel_utils.o:  $(BESSEL_UTILS_DIR)/bessel_utils.f90 bessel_functions_toms.o
	$(FC) $(FLAGS) -c  $(BESSEL_UTILS_DIR)/bessel_utils.f90

input_param_utils.o:  $(INPUT_PARAMS_DIR)/input_param_utils.f90 
	$(FC) $(FLAGS) -c  $(INPUT_PARAMS_DIR)/input_param_utils.f90

io_utils.o:  $(IO_UTILS_DIR)/io_utils.f90 input_param_utils.o
	$(FC) $(FLAGS) -c  $(IO_UTILS_DIR)/io_utils.f90

contour_init_utils.o: $(CONTOUR_UTILS_DIR)/contour_init_utils.f90 input_param_utils.o
	$(FC) $(FLAGS) -c  $(CONTOUR_UTILS_DIR)/contour_init_utils.f90

contour_generate_utils.o: $(CONTOUR_UTILS_DIR)/contour_generate_utils.f90 contour_init_utils.o input_param_utils.o
	$(FC) $(FLAGS) -c  $(CONTOUR_UTILS_DIR)/contour_generate_utils.f90

main.o: main.f90 bessel_utils.o io_utils.o input_param_utils.o contour_init_utils.o contour_generate_utils.o
	$(FC) $(FLAGS) -c main.f90

clean:
	rm -f *.o *.mod $(TARGET)

