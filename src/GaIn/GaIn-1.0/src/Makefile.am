lib_LIBRARIES = libGaIn.a

CFLAGS += -Wall

libGaIn_a_SOURCES = \
	CenteredPowerExpansion.f90 Coulomb.f90 \
	CoulombUtil.f90 CoulombUtilExpert.f90 CubicHarmonicProduct.f90 \
	Laplacian.f90 Overlap.f90 \
	R_1_norm.f90 R_Laplacian.f90 \
	R_from_Y.f90 R_to_D_Conversion.f90 \
	XYZ_power_to_ir.f90 Y_Value.f90 \
	Electrostatics.f90 R2Moment.f90 \
	Transform.f90 A_functions.h rlYlmLaplacian.f90

#temporary compiling line for gfortran
AM_FCFLAGS=-I./specialized_D_Coulomb_X -I$(top_srcdir)
AM_CFLAGS=-I$(top_builddir)

GaIn_HEADERS= GaIn.h GaIn.config.h
GaIndir= $(prefix)/include
CLEANFILES = *.mod *.MOD

#dependencies
Transform.o : Y_Value.o
CubicHarmonicProduct.o : CenteredPowerExpansion.o
Laplacian.o : R_1_norm.o R_Laplacian.o XYZ_power_to_ir.o R_from_Y.o CubicHarmonicProduct.o rlYlmLaplacian.o
Overlap.o : R_1_norm.o XYZ_power_to_ir.o R_from_Y.o CubicHarmonicProduct.o 
R2Moment.o : R_1_norm.o XYZ_power_to_ir.o R_from_Y.o CubicHarmonicProduct.o 
Coulomb.o : CoulombUtil.o CoulombUtilExpert.o R_to_D_Conversion.o XYZ_power_to_ir.o R_from_Y.o R_1_norm.o CubicHarmonicProduct.o
Electrostatics.o : CoulombUtil.o R_to_D_Conversion.o XYZ_power_to_ir.o R_from_Y.o R_1_norm.o CubicHarmonicProduct.o

SUBDIRS = specialized_D_Coulomb_X

libGaIn.a: $(libGaIn_a_OBJECTS) specialized_D_Coulomb_X/*.o
	$(AR) $(ARFLAGS) $@ $(libGaIn_a_OBJECTS) specialized_D_Coulomb_X/*.o
	$(RANLIB) $@

