.SUFFIXES:
.SUFFIXES:.o .f90 .f .mod

SHELL=/bin/sh

#for gfortran
#F90      =  gfortran -g
#F90FLAGS = -cpp -fbounds-check -fdefault-real-8 -O3  
#F90LINK  = $(F90FLAGS) -L. -lnspcg -L /DATA/ns_2D_3D_new_grid/CODE_16/CODE_NS_VISCOELASTIC_5/CODE_14 -O3


F90      =  gfortran -g
F90FLAGS = -cpp -fbounds-check  -fdefault-real-8 
F90LINK  = $(F90FLAGS) -L. liblapack.a librefblas.a -lnspcg -L ~/Documents/xzheng/DNS/CODE_14 -O3


##for ifort
#F90      =  ifort
#F90FLAGS = -fpp  #-O3   
#F90LINK  = $(F90FLAGS) -L. liblapack.a librefblas.a -lnspcg -L ~/Documents/xzheng/DNS/CODE_14 -O3

SOURCES = mBase.f90                \
	  mVitesse.f90             \
	  mOutils.f90              \
	  mPression.f90            \
	  mSolver.f90              \
	  mTimestep.f90            \
	  mChampsprof.f90          \
	  mPoste_traitement.f90    \
	  mHyperbolic_part.f90     \
	  mViscoelastic_models.f90 \

CODE_09_E1_BDF2_ABCN_RK3: $(SOURCES:.f90=.o) CODE_09_E1_BDF2_ABCN_RK3.f90
	$(F90) $(SOURCES:.f90=.o) CODE_09_E1_BDF2_ABCN_RK3.f90 $(F90LINK) -o CODE_09_E1_BDF2_ABCN_RK3.out
	./CODE_09_E1_BDF2_ABCN_RK3.out 
 

.f90.o :
	$(F90) $(F90FLAGS) -c $*.f90 -o $*.o

clean :
	
	\rm -f *.il *.mod *.o *.out *~ 
dep .depend:
	makedepf90 $(SOURCES) > .depend

include .depend





 
