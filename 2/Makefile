MPIF90 = mpif90 -Wall -Wextra -g -O3
GFORTRAN = gfortran -Wall -Wextra -g -O3

# the objects for test1

MPIOBJ = ModCommunication.o ModAdvance.o ModSavePlot.o
GFORTRANOBJ = 	ModYinYang.o ModConst.o ModAllocation.o	ModMath.o 	ModDeviation.o 	ModParameter.o \
		ModControl.o	ModSSM.o		ModSpherical.o	ModBlock.o \
		ModYinYangTree.o   	ModBoundary.o	ModLinReconstruct.o ModDiffusion.o ModEquation.o	ModGC.o \
		ModTimeStep.o

# # the objects for test2

TEST2OBJ = ModLinReconstruct.o

OBJ = $(GFORTRANOBJ) $(MPIOBJ)

MPISRC = ModCommunication.f90 ModAdvance.f90 ModSavePlot.f90
GFORTRANSRC = ModAllocation.f90 ModBlock.f90 ModBoundary.f90 ModDeviation.f90 \
              ModEquation.f90 ModGC.f90 ModMath.f90 ModOcTree.f90 ModParameter.f90 \
              ModSSM.f90 ModTimeStep.f90 ModVariable.f90
PROGRAM = test1.f90
TEST2_PROGRAM = test2.f90

EXE = test1.exe
EXE_test2 = test2.exe

%.o: %.f90
	$(MPIF90) -c $< -o $@

$(MPIOBJ): %.o: %.f90
	$(MPIF90) -c $< -o $@

#$(TEST2OBJ):  %.o: %.f90
#	$(GFORTRAN) -c $< -o $@

$(GFORTRANOBJ): %.o: %.f90
	$(GFORTRAN) -c $< -o $@

$(EXE): $(GFORTRANOBJ) $(MPIOBJ) $(PROGRAM)
	$(MPIF90) -o $(EXE) $(OBJ) $(PROGRAM)

$(EXE_test2): $(TEST2OBJ) $(TEST2_PROGRAM)
	$(GFORTRAN) -o $(EXE_test2) $(TEST2OBJ) $(TEST2_PROGRAM)

clean:
	rm -f $(OBJ) $(EXE)
