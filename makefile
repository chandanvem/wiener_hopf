LIBNAG           =  /usr/local/lib32/NAG/fllux23dcl/lib/libnag_nag.a
FORT             =  ifort
OPTF             = -O3 -Dintel_ -DALLOW_NON_INIT -nofor_main
OPTL             = -O3 -nofor_main
OPTS             = -O3
#FORT             =  nagfor
#OPTF             = -O4 -fpp -Dintel_ -DALLOW_NON_INIT -kind=byte
#OPTL             = -O4
DEBUG            = 

OBJECTS =  main.f90

main: 
	$(FORT) $(OPTS) $(OBJECTS) $(LIBNAG) -o main

clean:
	rm *.o -f
	rm le -f
	rm -rf main
	rm -rf *.out
	rm -rf mesh_polar.*
