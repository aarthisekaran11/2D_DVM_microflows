AUTOMAKE_OPTIONS = foreign
SUBDIRS          = src test

# List all the input files we want installed
input_file_SRC = $(top_srcdir)/input_files/DVM_grvy.inp

# Make the input_files installed in $(prefix)/input_files directory
input_filedir = $(prefix)/input_files
input_file_DATA  = $(input_file_SRC)

# This is what gets distributed when we "make dist"
EXTRA_DIST       = AUTHORS CHANGES README COPYING LICENSE \
                   $(input_file_SRC)

AUX_DIST         = build-aux/install-sh
AUX_DIST        += build-aux/missing
AUX_DIST        += build-aux/config.guess
AUX_DIST        += build-aux/config.sub
AUX_DIST        += build-aux/depcomp
AUX_DIST        += build-aux/ltmain.sh

# Files to be deleted by 'make maintainer-clean'
MAINTAINERCLEANFILES = aclocal.m4                               \
                       aminclude.am                             \
                       autom4te.cache/*                         \
                       $(AUX_DIST)                              \
                       config.log                               \
                       config.status                            \
                       config.sub                               \
                       configure                                \
                       config.h                                 \
                       config.h.in                              \
                       Makefile.in                              \
                       src/Makefile.in                          \
                       test/Makefile.in

.PHONY: model_doc
model_doc:
	@cd docs/modeling_document; make

.PHONY: clean_model_doc
clean_model_doc:
	@cd docs/modeling_document; make clean
