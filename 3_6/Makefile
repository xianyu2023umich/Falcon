MPIF90		=	mpif90 -Wall -Wextra -g -O3
GFORTRAN 	= 	gfortran -Wall -Wextra -g -O3 

test0_PGM	=	test0.exe
test0_OBJ   =   ModMath.o

test1_PGM	=	test1.exe
test1_MPIOBJ=	ModCommunication.o 	ModCheck.o \
				ModAdvance.o		
test1_OBJ   =   ModYinYang.o		ModDerivative.o		ModSpherical.o		 \
				ModConst.o			ModParameters.o		ModLinReconstruct.o ModLookUpTable.o  \
				ModEOS.o			ModOpacity.o		ModMath.o 			ModStratification_new.o\
				ModAllocation.o	    ModControl.o		ModBlock.o \
				ModBoundary.o		ModDiffusion.o      ModDivB.o			ModEquation.o       	\
				ModYinYangTree.o   	ModAMR.o			ModInitiation.o		\
			    ModReadParameters.o	ModGC.o			    ModWaveSpeed.o		ModTimeStep.o		ModSavePlot.o

test_CH_PGM =   test_coronalheating.exe
test_CH_OBJ =   ModConst.o          ModParameters.o     ModStratification.o ModBlock.o          \
                ModMHD.o           ModDerivative.o     ModSpherical.o      ModCoronalHeating.o 

test_lookuptable_PGM =   test_lookuptable.exe
test_lookuptable_OBJ =   ModLookUpTable.o

test_opacity_PGM =   test_opacity.exe
test_opacity_OBJ =   ModLookUpTable.o ModOpacity.o

test_eos_PGM =   test_eos.exe
test_eos_OBJ =   ModLookUpTable.o ModEOS.o

test_stratification_new_PGM =   test_stratification_new.exe
test_stratification_new_OBJ =   ModLookUpTable.o ModConst.o ModMath.o ModEOS.o ModOpacity.o ModStratification_new.o

test_stratification_new_real_PGM =   test_stratification_new_real.exe
test_stratification_new_real_OBJ =   ModLookUpTable.o ModConst.o ModMath.o ModEOS.o ModOpacity.o ModStratification_new.o

test_ModEquation_PGM =   test_ModEquation.exe
test_ModEquation_OBJ =   ModConst.o ModLookUpTable.o ModMath.o ModEOS.o ModOpacity.o ModParameters.o \
					ModStratification_new.o ModBlock.o ModBoundary.o ModDerivative.o ModSpherical.o \
					ModLinReconstruct.o ModDiffusion.o ModEquation.o

test_ModAdvance_PGM =   test_ModAdvance.exe
test_ModAdvance_OBJ =   ModConst.o ModLookUpTable.o ModMath.o ModEOS.o ModOpacity.o ModParameters.o \
				ModStratification.o ModStratification_new.o ModControl.o ModAllocation.o ModYinYang.o \
				ModBlock.o ModYinYangTree.o ModGC.o ModBoundary.o ModDerivative.o ModSpherical.o ModLinReconstruct.o \
				ModDiffusion.o ModEquation.o ModTimeStep.o ModCommunication.o ModCheck.o \
				ModReadParameters.o ModAdvance.o

test_ModSavePlot_PGM =  test_ModSavePlot.exe
test_ModSavePlot_OBJ =   ModConst.o ModLookUpTable.o ModMath.o ModEOS.o ModOpacity.o ModParameters.o \
				ModStratification.o ModStratification_new.o ModControl.o ModAllocation.o ModYinYang.o \
				ModBlock.o ModYinYangTree.o ModGC.o ModBoundary.o ModDerivative.o ModSpherical.o ModLinReconstruct.o \
				ModDiffusion.o ModEquation.o ModTimeStep.o ModCommunication.o ModCheck.o \
				ModReadParameters.o ModAMR.o ModInitiation.o ModSavePlot.o

test_ModAMR_PGM =   test_ModAMR.exe
test_ModAMR_OBJ =   ModConst.o ModLookUpTable.o ModMath.o ModEOS.o ModOpacity.o ModParameters.o \
				ModStratification.o ModStratification_new.o ModControl.o ModAllocation.o ModYinYang.o \
				ModBlock.o ModYinYangTree.o ModGC.o ModBoundary.o ModDerivative.o ModSpherical.o ModLinReconstruct.o \
				ModDiffusion.o ModEquation.o ModTimeStep.o ModCommunication.o ModCheck.o \
				ModReadParameters.o ModAMR.o

test_ModInitiation_PGM =   test_ModInitiation.exe
test_ModInitiation_OBJ =   ModConst.o ModLookUpTable.o ModMath.o ModEOS.o ModOpacity.o ModParameters.o \
				ModStratification.o ModStratification_new.o ModControl.o ModAllocation.o ModYinYang.o \
				ModBlock.o ModYinYangTree.o ModGC.o ModBoundary.o ModDerivative.o ModSpherical.o ModLinReconstruct.o \
				ModDiffusion.o ModEquation.o ModTimeStep.o ModCommunication.o ModCheck.o \
				ModReadParameters.o ModAMR.o ModInitiation.o ModSavePlot.o


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

test_CH: $(test_CH_OBJ) test_coronalheating.f90
	$(GFORTRAN)	-o	$(test_CH_PGM)	$(test_CH_OBJ) test_coronalheating.f90

#$(test_CH_OBJ)	:	%.o	:	%.f90
#	$(GFORTRAN)	-c 	$< 	-o 	$@

test_lookuptable: $(test_lookuptable_OBJ) test_lookuptable.f90
	$(GFORTRAN)	-o	$(test_lookuptable_PGM)	$(test_lookuptable_OBJ) test_lookuptable.f90

test_opacity: $(test_opacity_OBJ) test_opacity.f90
	$(GFORTRAN)	-o	$(test_opacity_PGM)	$(test_opacity_OBJ) test_opacity.f90

test_eos: $(test_eos_OBJ) test_eos.f90
	$(GFORTRAN)	-o	$(test_eos_PGM)	$(test_eos_OBJ) test_eos.f90

test_stratification_new: $(test_stratification_new_OBJ) test_stratification_new.f90
	$(GFORTRAN)	-o	$(test_stratification_new_PGM)	$(test_stratification_new_OBJ) test_stratification_new.f90

test_stratification_new_real: $(test_stratification_new_real_OBJ) test_stratification_new_real.f90
	$(GFORTRAN)	-o	$(test_stratification_new_real_PGM)	$(test_stratification_new_real_OBJ) test_stratification_new_real.f90

test_ModEquation: $(test_ModEquation_OBJ) test_ModEquation.f90
	$(GFORTRAN)	-o	$(test_ModEquation_PGM)	$(test_ModEquation_OBJ) test_ModEquation.f90

test_ModAdvance: $(test_ModAdvance_OBJ) test_ModAdvance.f90
	$(MPIF90)	-o	$(test_ModAdvance_PGM)	$(test_ModAdvance_OBJ) test_ModAdvance.f90

test_ModSavePlot: $(test_ModSavePlot_OBJ) test_ModSavePlot.f90
	$(MPIF90)	-o	$(test_ModSavePlot_PGM)	$(test_ModSavePlot_OBJ) test_ModSavePlot.f90

test_ModAMR: $(test_ModAMR_OBJ) test_ModAMR.f90
	$(MPIF90)	-o	$(test_ModAMR_PGM)	$(test_ModAMR_OBJ) test_ModAMR.f90

test_ModInitiation: $(test_ModInitiation_OBJ) test_ModInitiation.f90
	$(MPIF90)	-o	$(test_ModInitiation_PGM)	$(test_ModInitiation_OBJ) test_ModInitiation.f90

test_ModCommunication_compile: ModConst.o ModLookUpTable.o ModMath.o ModEOS.o ModOpacity.o ModParameters.o \
				ModStratification.o ModStratification_new.o ModControl.o ModAllocation.o ModYinYang.o \
				ModBlock.o ModYinYangTree.o ModGC.o ModCommunication.o
	@echo "ModCommunication compiled successfully."

ModSavePlot.o: ModSavePlot.f90
	$(MPIF90)	-c 	$< 	-o 	$@

%.o: %.f90
	$(GFORTRAN)	-c 	$< 	-o 	$@

test_divB_cleaning: test_divB_cleaning.f90
	$(GFORTRAN) -o test_divB_cleaning.exe test_divB_cleaning.f90

test_divB_sph: test_divB_sph.f90
	$(GFORTRAN) -o test_divB_sph.exe test_divB_sph.f90

test_divB_sph_validate: test_divB_sph_validate.f90
	$(GFORTRAN) -o test_divB_sph_validate.exe test_divB_sph_validate.f90

clean:
	rm -f *.o *.exe *.mod
