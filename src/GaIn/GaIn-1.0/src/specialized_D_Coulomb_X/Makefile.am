noinst_LIBRARIES = lib_spec_D_Coulomb_X.a

lib_spec_D_Coulomb_X_a_SOURCES = \
	my_erf.c \
	A_functions.f90 \
	Av_functions.h \
	D_X_R_0.f90 \
	D_X_R_1.f90 \
	D_X_R_2.f90 \
	D_X_R_3.f90 \
	D_X_R_4.f90 \
	D_X_Y_0.f90 \
	D_X_Y_1.f90 \
	D_X_Y_2.f90 \
	D_X_Y_3.f90 \
	D_X_Y_4.f90 \
	D_X_R_shell_0.f90 \
	D_X_R_shell_1.f90 \
	D_X_R_shell_2.f90 \
	D_X_R_shell_3.f90 \
	D_X_R_shell_4.f90 \
	D_X_Y_shell_0.f90 \
	D_X_Y_shell_1.f90 \
	D_X_Y_shell_2.f90 \
	D_X_Y_shell_3.f90 \
	D_X_Y_shell_4.f90 \
	D_X_R_shell_0_vec1.f90 \
	D_X_R_shell_1_vec1.f90 \
	D_X_R_shell_2_vec1.f90 \
	D_X_R_shell_3_vec1.f90 \
	D_X_R_shell_4_vec1.f90 \
	D_X_Y_shell_0_vec1.f90 \
	D_X_Y_shell_1_vec1.f90 \
	D_X_Y_shell_2_vec1.f90 \
	D_X_Y_shell_3_vec1.f90 \
	D_X_Y_shell_4_vec1.f90 \
	D_X_R_shell_0_vec4.f90 \
	D_X_R_shell_1_vec4.f90 \
	D_X_R_shell_2_vec4.f90 \
	D_X_R_shell_3_vec4.f90 \
	D_X_R_shell_4_vec4.f90 \
	D_X_Y_shell_0_vec4.f90 \
	D_X_Y_shell_1_vec4.f90 \
	D_X_Y_shell_2_vec4.f90 \
	D_X_Y_shell_3_vec4.f90 \
	D_X_Y_shell_4_vec4.f90 \
	D_X_V_and_E.f90 \
	D_X_V_and_E_vec1.f90 \
	D_X_V_and_E_vec4.f90

#temporary compiling line for gfortran
AM_FCFLAGS = 
CLEANFILES = *.mod *.MOD

#dependencies
my_erf.o : my_erf.c
A_functions.o : A_functions.f90
$(lib_spec_D_Coulomb_X_a_OBJECTS): $(lib_spec_D_Coulomb_X_a_SOURCES) A_functions.o
lib_spec_D_Coulomb_X.a: $(lib_spec_D_Coulomb_X_a_OBJECTS)
	touch $@
