
 PROG = ~/bin/bff

 SOURCES=\
 modules.f90\
 version.f90\
 eval_opt.f90\
 io.f90\
 string.f90\
 math.f90\
 LennardJones.f90\
 bond.f90\
 param.f90\
 print.f90\
 atype.f90\
 numhess.f90\
 bonded.f90\
 printmat.f90\

# sources needed pre-processing
 FSOURCES=\
 main.F90\

OBJS=$(SOURCES:.f90=.o) $(FSOURCES:.F90=.o)

BUILID:=$(shell date)
GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always --tags)
$(info Building: $(GIT_VERSION))


#------- openmp parallelization --------------
# set to "yes" if openmp is desired
# "make clean" is needed after (de)activation!
USE_OMP=no
#------- openmp parallelization --------------

OPENBLAS=/usr/qc/openblas_lib/


# ***   GFORTRAN ***
  FC = gfortran -static
  FLAGS= -O -ffree-line-length-none -m64
# FLAGS= -Og -g -fbounds-check -ffree-line-length-none -m64
  LIBS= -llapack -lblas
#  LIBS= -L$(OPENBLAS)/lib/ -lopenblas -lpthread


# ***   PGI       ***
# FC=pgfortran -Bstatic
# FLAGS= -fast -Minfo=ipa,opt,par,mp,inline -i8 -traceback
# FLAGS= -fastsse -Mipa=fast  -Minfo=opt,par,loop,inline,vect,mp -i8 -traceback -Mconcur=levels:4 -Mvect=levels:4
# LIBS= -llapack -lblas


ifeq ($(USE_OMP),yes)
 DFLAGS= -DOPENMP
# FLAGS+= -mp
 FC+= -fopenmp
endif

# targets:
.PHONY: all
.PHONY: clean
.PHONY: archive

all: version $(PROG)


version:
	@touch version.f90
	@echo 'writing new version.f90'
	@echo "subroutine version" > version.f90
	@echo     " print*,  'Build info:'" >> version.f90
#	@echo     " print*,  ' github        : https://github.com/hokru/struca '">> version.f90
	@echo     " print*,  ' build date    : $(BUILID)     '"      >> version.f90
	@echo     " print*,  ' git version   : $(GIT_VERSION)'"      >> version.f90
	@echo     " print*,  ' '          "     >> version.f90
	@echo "end subroutine" >> version.f90


%.o: %.f90
	@echo "making $@ from $<"
	$(FC) $(FLAGS) -c $< -o $@

%.o: %.F90
	@echo "making $@ from $<"
	$(FC) $(DFLAGS) $(FLAGS)  -c $< -o $@

$(PROG):$(OBJS)
	$(FC) $(LINK) $(OBJS) $(LIBS) -o $(PROG)


clean:
	rm -f *.o *.mod $(PROG)

archive:
	git archive -o archive_$(GIT_VERSION).zip HEAD
