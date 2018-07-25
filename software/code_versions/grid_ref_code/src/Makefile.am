BUILT_SOURCES   = .license.stamp

AUTOMAKE_OPTIONS = foreign
SUBDIRS = .

all-local: Makefile.dep1 Makefile.dep2

#----------------------------------------
# Programs and libraries we want to build
#----------------------------------------
# Idea here is to have dvm be the executable that's
# built from main.f90 and the precompiled library
# FIXME: Change this name to something better
bin_PROGRAMS = dvm_exe

# This is the library that will have everything else
# FIXME: Change this name to something better

# Karl suggested that the static libraries will get correct parallel build
#lib_LTLIBRARIES = libdvm.la
noinst_LIBRARIES = libdvm.a

#----------------------------------------
# List all the source files here
#----------------------------------------
cfiles_SOURCES = $(top_srcdir)/src/c_files/get_gsl_rng_default.c \
                 $(top_srcdir)/src/c_files/polyfit.c             \
                 $(top_srcdir)/src/c_files/lnfit.c               \
                 $(top_srcdir)/src/c_files/temp_vib_fit.c        \
                 $(top_srcdir)/src/c_files/rotfit.c              \
                 $(top_srcdir)/src/c_files/temp_rot_fit.c

utils_SOURCES = $(top_srcdir)/src/utils/constants.f90                \
                $(top_srcdir)/src/utils/error_check.f90              \
                $(top_srcdir)/src/utils/conversion.f90               \
                $(top_srcdir)/src/utils/random_number_generation.f90 \
                $(top_srcdir)/src/utils/math_utils.f90               \
                $(top_srcdir)/src/utils/sort_and_search.f90


data_SOURCES = $(top_srcdir)/src/input_output/species_data.f90

scaling_SOURCES = $(top_srcdir)/src/solver/time_stepping.f90 \
                  $(top_srcdir)/src/solver/scaling.f90

grid_SOURCES = $(top_srcdir)/src/grid/velocity_grid.f90 \
	       $(top_srcdir)/src/grid/physical_grid.f90

dist_funcs_SOURCES = $(top_srcdir)/src/dist_funcs/dist_func.f90      \
		     $(top_srcdir)/src/dist_funcs/physical_props.f90

conditions_SOURCES = $(top_srcdir)/src/conditions/boundary_conditions.f90 \
                     $(top_srcdir)/src/conditions/shock_conditions.f90    \
		     $(top_srcdir)/src/conditions/initial_conditions.f90 

collisions_SOURCES = $(top_srcdir)/src/collisions/collision_utilities.f90     \
                     $(top_srcdir)/src/collisions/replenishing_collisions.f90 \
                     $(top_srcdir)/src/collisions/full_mc.f90                 \
                     $(top_srcdir)/src/collisions/variance_reduction.f90      \
                     $(top_srcdir)/src/collisions/collision_integral.f90

convection_SOURCES = $(top_srcdir)/src/convection/convection.f90

input_output_SOURCES = $(top_srcdir)/src/input_output/restart.f90       \
                       $(top_srcdir)/src/input_output/screen_output.f90 \
                       $(top_srcdir)/src/input_output/command_line.f90  \
                       $(top_srcdir)/src/input_output/visualization.f90 \
                       $(top_srcdir)/src/input_output/read_input.f90

solver_SOURCES = $(top_srcdir)/src/solver/dvm_solver.f90


dvm_SRC = $(cfiles_SOURCES)       \
          $(utils_SOURCES)        \
          $(grid_SOURCES)         \
          $(scaling_SOURCES)      \
          $(data_SOURCES)         \
          $(dist_funcs_SOURCES)   \
          $(top_srcdir)/src/collisions/remapping.f90 \
	  $(conditions_SOURCES)   \
          $(collisions_SOURCES)   \
          $(convection_SOURCES)   \
          $(input_output_SOURCES) \
          $(solver_SOURCES)

libdvm_a_SOURCES =  $(dvm_SRC)

dvm_exe_SOURCES = $(top_srcdir)/src/main.f90

#AM_FCFLAGS = -I$(GRVY_PREFIX)/lib
#dvm_exe_LDFLAGS = -L$(top_builddir)/src -L$(GRVY_PREFIX)/lib
#dvm_exe_LDADD   = -ldvm -lgrvy #-lgrvyf

# switch to using grvy macro variables
AM_CPPFLAGS = $(GSL_CFLAGS)
AM_LDFLAGS = 
AM_FCFLAGS    = $(GRVY_FCFLAGS)
dvm_exe_LDADD = -ldvm $(GRVY_FLIBS) $(GSL_LIBS) -L$(top_builddir)/src


# ----------------------------------------
# Automatic dependency generation for 
# Fortran (libtool based dependencies)
# ----------------------------------------
IGNORED_MODULES = -i grvy -i iso_c_binding
PARSE_OPTIONS = -I $(GRVY_PREFIX)/lib -I $(top_srcdir) $(IGNORED_MODULES) -strip

Makefile.dep1: $(libdvm_a_SOURCES)
	$(top_srcdir)/build-aux/common/sfmakedepend.pl $(PARSE_OPTIONS) $^ > $@

Makefile.dep2: $(dvm_exe_SOURCES)
	$(top_srcdir)/build-aux/common/sfmakedepend.pl $(PARSE_OPTIONS) $^ > $@

#This forces the dependencies generated to be included
#in Automakes dependency graph
-include Makefile.dep1
-include Makefile.dep2

#----------------------------
# Extra data for installation
#----------------------------
moduledir = $(prefix)/lib
module_DATA =         

#Files to be deleted by 'make clean'
CLEANFILES = *.mod Makefile.dep1 Makefile.dep2     \
             $(top_srcdir)/src/collisions/*.o      \
             $(top_srcdir)/src/collisions/*.lo     \
             $(top_srcdir)/src/dist_funcs/*.o      \
             $(top_srcdir)/src/dist_funcs/*.lo     \
             $(top_srcdir)/src/input_output/*.o    \
             $(top_srcdir)/src/input_output/*.lo   \
             $(top_srcdir)/src/utils/*.o           \
             $(top_srcdir)/src/utils/*.lo

#---------------------------------
# Embedded license header support
#--------------------------------- 
.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/common/lic_utils/update_license.pl $(top_srcdir)/LICENSE $(dvm_SRC) $(dvm_exe_SOURCES)
	echo 'updated source license headers' >$@

