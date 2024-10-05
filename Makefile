MPIF90 = mpif90 -Wall -Wextra -g
GFORTRAN = gfortran -Wall -Wextra -g


MPIOBJ = ModCommunication.o ModAdvance.o ModSavePlot.o
GFORTRANOBJ = 	ModAllocation.o	ModMath.o 	ModDeviation.o 	ModParameter.o \
		ModControl.o	ModSSM.o	ModVariable.o	ModBlock.o \
		ModOcTree.o   	ModBoundary.o	ModEquation.o	ModGC.o \
		ModTimeStep.o

OBJ = $(GFORTRANOBJ) $(MPIOBJ)

MPISRC = ModCommunication.f90 ModAdvance.f90 ModSavePlot.f90
GFORTRANSRC = ModAllocation.f90 ModBlock.f90 ModBoundary.f90 ModDeviation.f90 \
              ModEquation.f90 ModGC.f90 ModMath.f90 ModOcTree.f90 ModParameter.f90 \
              ModSSM.f90 ModTimeStep.f90 ModVariable.f90
PROGRAM = test1.f90

EXE = test1.exe

%.o: %.f90
	$(MPIF90) -c $< -o $@

$(MPIOBJ): %.o: %.f90
	$(MPIF90) -c $< -o $@

$(GFORTRANOBJ): %.o: %.f90
	$(GFORTRAN) -c $< -o $@

$(EXE): $(GFORTRANOBJ) $(MPIOBJ) $(PROGRAM)
	$(MPIF90) -o $(EXE) $(OBJ) $(PROGRAM)

clean:
	rm -f $(OBJ) $(EXE)
