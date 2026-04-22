#
FC=gfortran

FV:	FV_ext.o FV_scal.F90	
	$(FC)	FV_scal.F90	FV_ext.o

FV_ext.o : FV_ext.F90
	$(FC) -c FV_ext.F90