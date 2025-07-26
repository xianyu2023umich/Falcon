MPIF90		=	mpif90 -Wall -Wextra -g -O3
GFORTRAN 	= 	gfortran -Wall -Wextra -g -O3

test0_PGM	=	test0.exe
test0_OBJ   =   ModMath.o

test1_PGM	=	test1.exe
test1_MPIOBJ=	ModCommunication.o 	ModLinearSolver.o	ModPFSS.o			ModAdvance.o \
				ModSavePlot.o		ModCheck.o
test1_OBJ   =   ModVariables.o 		ModYinYang.o		ModDerivative.o		ModSpherical.o		 \
				ModConst.o			ModParameters.o		ModLinReconstruct.o ModStratification.o \
				ModMath.o			ModAllocation.o	    ModControl.o		ModBlock.o \
				ModBoundary.o		ModDiffusion.o      ModEquation.o       	\
				ModYinYangTree.o   	ModAMR.o			ModInitiation.o		ModMultigrid.o		\
			    ModReadParameters.o	ModGC.o			    ModTimeStep.o		ModMagnetogram.o



#$(test0_OBJ)	:	%.o	:	%.f90
#	$(GFORTRAN)	-c 	$< 	-o 	$@

test0: $(test0_OBJ) test0.f90
	$(GFORTRAN)	-o	$(test0_PGM)	$(test0_OBJ) test0.f90

$(test1_OBJ)	:	%.o	:	%.f90
	$(GFORTRAN)	-c 	$< 	-o 	$@

$(test1_MPIOBJ)	:	%.o	:	%.f90
	$(MPIF90)	-c 	$< 	-o 	$@

test1: $(test1_OBJ) $(test1_MPIOBJ) test1.f90
	$(MPIF90)	-o	$(test1_PGM)	$(test1_OBJ) $(test1_MPIOBJ) test1.f90

clean:
	rm -f *.o *.exe *.mod