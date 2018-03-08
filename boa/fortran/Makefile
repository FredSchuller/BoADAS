verbosity = no

# Machine specific variables are set in file Makevars:
include Makevars

# standard modules
modules = fFlag fStat fFlagHandler fWavelets

# modules which need a special flag for non optimization
mod_opt = fUtilities

# module which need lapack
mod_lapack = fFit 

# module which needs lapack + non-opt
mod_opt_lapack = fMap

# module which need fit
mod_fit = fSNF fBaseline

#----------------------------------------------------------

all: $(modules) $(mod_opt) $(mod_svd) $(mod_lapack) $(mod_fit) $(mod_opt_lapack)
#	chcon -t texrel_shlib_t *.so


FF            = gfortran
F2PY_compiler = --f90exec=gfortran --f77exec=gfortran
F2PY_extra    = 

ifeq "$(verbosity)" "yes"
  F2PY_verbose = --verbose -DF2PY_REPORT_ON_ARRAY_COPY=50000 
else
  F2PY_verbose = --quiet
endif

F2PY        = f2py -c
F2PY_option = $(F2PY_compiler) --build-dir . -m

.SUFFIXES:
.SUFFIXES: .f90 .o .so .mod


$(modules): 
	$(F2PY) $(F2PY_verbose) $(F2PY_option) $@ $@.f90 $(F2PY_libs)
$(mod_opt):
	$(F2PY) $(F2PY_verbose) $(F2PY_extra) $(F2PY_option) $@ $@.f90 $(F2PY_libs)
$(mod_fit): fFit.o
	$(F2PY) $(F2PY_verbose) $(F2PY_option) --link-lapack_opt $@ $@.f90 $? $(F2PY_libs)
$(mod_lapack): 
	$(F2PY) $(F2PY_verbose) $(F2PY_option) --link-lapack_opt $@ $@.f90 $? $(F2PY_libs)
$(mod_opt_lapack):
	$(F2PY) $(F2PY_verbose) $(F2PY_extra) $(F2PY_option) --link-lapack_opt $@ $@.f90 $? $(F2PY_libs)

.PHONY : clean

clean:
	rm -f *.o *.mod
	rm -rf src

cleanall: clean
	rm -f *.so
