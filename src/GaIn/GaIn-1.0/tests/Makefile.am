check_PROGRAMS = test_transform test_overlap test_r2moment test_laplacian test_coulomb test_electrostatics $(extra_tests)
extra_tests = test_ir_index test_R_1_norm

LIBS += -lGaIn
AM_LDFLAGS= -L$(top_builddir)/src/
AM_FCFLAGS= -I$(top_builddir)/src/
test_coulomb_SOURCES= test_Coulomb.f90
test_electrostatics_SOURCES=test_Electrostatics.f90
test_overlap_SOURCES= test_Overlap.f90
test_r2moment_SOURCES= test_R2Moment.f90
test_laplacian_SOURCES= test_Laplacian.f90
test_transform_SOURCES= test_Transform.f90

test_ir_index_SOURCES= test_ir_index.f90
test_R_1_norm_SOURCES= test_R_1_norm.f90

check: 
	./test_transform 
	./test_ir_index 
	./test_R_1_norm 
	./test_overlap
	./test_r2moment
	./test_laplacian
	./test_electrostatics
	./test_coulomb
