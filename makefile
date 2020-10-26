# Generate PTX use ‑Mcuda=keepptx

GCC=pgcc
FC=pgfortran
#GCFLAGS_ACC=-O4 -acc -Mpreprocess -Minfo=accel -mp=allcores -Minline -DITHREADS=16 -DJTHREADS=8 -ta=tesla:keepptx -tp=penryn
#GCFLAGS_ACC=-O4 -acc -Mpreprocess -Minfo=accel -mp=allcores -Minline -DITHREADS=16 -DJTHREADS=8 -ta=radeon:keep -tp=penryn
#GCFLAGS_ACC=-O4 -acc -Mpreprocess -Minfo=accel -mp=allcores -Minline -DITHREADS=16 -DJTHREADS=8 -ta=radeon:tahiti -tp=penryn
#GCFLAGS_ACC=-O4 -acc -Mpreprocess -Minfo=accel -mp=allcores -Minline -DITHREADS=16 -DJTHREADS=8 -ta=host -tp=penryn
#GCFLAGS_ACC=-O4 -acc -Mpreprocess -Minfo=accel -mp=allcores -Minline -DITHREADS=16 -DJTHREADS=8 -ta=tesla -tp=penryn
#GCFLAGS_ACC=-O4 -acc -Mpreprocess -Minfo=accel  -mp=allcores -Minline -DITHREADS=16 -DJTHREADS=8 -ta=tesla -tp=penryn -Mcuda=ptxinfo,Xptxas,v,abi=no -pgf90libs 
GCFLAGS_ACC=-O4 -acc -Mpreprocess -Minfo=accel -mp=allcores -Minline -DITHREADS=8 -DJTHREADS=16 -tp=penryn -Mcuda=ptxinfo -pgf90libs 
#-ta=tesla:cc35,cuda5.5,maxregcount:76  -Mcuda=ptxinfo -pgf90libs
#-ta=tesla:cc35,cuda5.5
#-Mcuda=ptxinfo -pgf90libs -ta=tesla:cc30,cuda5.5,maxregcount:76
#-ta=tesla:cc35,cuda5.5
#-ta=tesla:cc35,cuda6.0,maxregcount:76
#GCFLAGS_ACC=-O4 -acc -ta=host -ta=tesla:cc35 -Mpreprocess -Minfo=accel -mp=allcores -Minline -DITHREADS=16 -DJTHREADS=8 -tp=penryn
#GCFLAGS_ACC=-acc -ta=host -ta=tesla,cc35 -Mpreprocess -Minfo=accel -mp=allcores -Minline
#FCFLAGS=-O4 -Minfo -Mpreprocess -mp=all -Mvect=simd:256 -Minline #-traceback -openmp -parallel -fpp
FCFLAGS_ACC=-O4 -acc -ta=host  -Mpreprocess -Minfo=accel -mp=allcores -Minline  -DITHREADS=16 -DJTHREADS=8 -tp=penryn
FCFLAGS=-O4 -Mpreprocess -Minfo=accel -mp=allcores -Minline #-traceback -openmp -parallel -fpp
LDC_INCLUDE=-I functions/


all: ldc

#ldc: interface_f_2_c_vars.o set_precision.o setup.o matrix_manip.o solvers.o ldc.o
#	$(FC) $(LDC_INCLUDE) $(FCFLAGS) interface_f_2_c_vars.o set_precision.o setup.o fileio.o solver_c.o matrix_manip.o solvers.o ldc.o -o $@ -lm

ldc: setup.o fileio.o solver_c.o ldc.o
	$(GCC) -I ./ $(GCFLAGS_ACC) -lgomp setup.o fileio.o solver_c.o ldc.o -o $@ -lm

#ldc_acc: interface_f_2_c_vars.o set_precision.o setup.o matrix_manip.o solvers_acc.o ldc.o
#	$(FC) $(LDC_INCLUDE) $(FCFLAGS_ACC) interface_f_2_c_vars.o set_precision.o setup.o fileio.o solver_c.o matrix_manip.o solvers_acc.o ldc.o -o $@ -lm

#interface_f_2_c_vars.o: interface_f_2_c_vars.f90
#	$(FC) $(LDC_INCLUDE) $(FCFLAGS) -c interface_f_2_c_vars.f90

#set_precision.o: set_precision.f90
#	$(FC) $(LDC_INCLUDE) $(FCFLAGS) -c set_precision.f90


setup.o: setup.c
	$(GCC) $(GCFLAGS_ACC) -c setup.c

fileio.o: fileio.c
	$(GCC) $(GCFLAGS_ACC) -c fileio.c

solver_c.o: solver_c.c
	$(GCC) $(GCFLAGS_ACC) -c solver_c.c

ldc.o: setup.o fileio.o solver_c.o ldc.c
	$(GCC) $(GCFLAGS_ACC) -c ldc.c

#matrix_manip.o: set_precision.o matrix_manip.f90
#	$(FC) $(LDC_INCLUDE) $(FCFLAGS) -c matrix_manip.f90

#solvers_acc.o: set_precision.o interface_f_2_c_vars.o matrix_manip.o solvers.f90
#	$(FC) $(LDC_INCLUDE) $(FCFLAGS_ACC) -c solvers.f90 -o $@

#solvers.o: set_precision.o interface_f_2_c_vars.o solver_c.o matrix_manip.o solvers.f90
#	$(FC) $(LDC_INCLUDE) $(FCFLAGS) -c solvers.f90

#ldc.o: fileio.o interface_f_2_c_vars.o solvers.o ldc.f90
#	$(FC) $(LDC_INCLUDE) $(FCFLAGS) -c ldc.f90

clean:
	rm -f *.o *.mod ldc
