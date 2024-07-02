Machine = local
FC	= gfortran 
FFLAGS  = -fopenmp -O3 -w -fcheck=all -g 

ifeq (${Machine},nias)
fftw3lib_dir    = /scratch/Saurabh/FFTW3/install/lib64     #/home/me23jok/ProjectX/FFTW3/lib
LDFLAGS         = -L${fftw3lib_dir} -lfftw3 -lfftw3_omp -lm

else ifeq (${Machine},ara)
fftw3lib_dir    = /home/me23jok/ProjectX/FFTW3/lib
LDFLAGS         = -L${fftw3lib_dir} -lfftw3 -lfftw3_omp -lm

else ifeq (${Machine},draco)
fftw3lib_dir    = /home/me23jok/fftw-3.3.10/lib
LDFLAGS         = -L${fftw3lib_dir} -lfftw3 -lfftw3_omp -lm

else ifeq (${Machine},local)
LDFLAGS         = -lfftw3 -lfftw3_omp -lm -llapack -lblas
endif


default: 	dynamics

dynamics:   	main.o   adiabatic.o\
	propagation.o	setpot.o\
	nuclear_wv.o\
	pulse.o
				
		${FC} *.o ${FFLAGS} ${LDFLAGS} -o dynamics
		
main.o:	 main.f90
	${FC} main.f90 ${FFLAGS} -c -o main.o 

setpot.o:	setpot.f90
	${FC} setpot.f90 ${FFLAGS} -c -o setpot.o 	
		
adiabatic.o:	adiabatic.f90
	${FC} adiabatic.f90 ${FFLAGS} -c -o adiabatic.o 
	
propagation.o:	propagation.f90
	${FC} propagation.f90 ${FFLAGS} -c -o propagation.o 
		
nuclear_wv.o:	nuclear_wv.f90
	${FC} nuclear_wv.f90 ${FFLAGS} -c -o nuclear_wv.o

pulse.o:	pulse.f90
	${FC} pulse.f90 ${FFLAGS} -c -o pulse.o

clean:
	rm -f *.o 
	rm -f *.mod 
	rm -f dynamics 

