
 PROG = ~/bin/imff

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
 main.f90\

OBJS=$(SOURCES:.f90=.o)

BUILID:=$(shell date)
GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always --tags)
$(info Building: $(GIT_VERSION))


OPENBLAS=/usr/qc/openblas_lib/

  FC = gfortran 
  FLAGS= -O3 -ffree-line-length-none -m64 
#  FLAGS= -Og -g -fbounds-check -ffree-line-length-none -m64 
#  LIBS= -llapack -lblas
#  LIBS= -L$(OPENBLAS)/lib/ -lopenblas -lpthread

# targets:
.PHONY: all
.PHONY: clean

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


$(PROG):$(OBJS) 
	$(FC) $(LINK) $(OBJS) $(LIBS) -o $(PROG)





clean:
	rm -f *.o *.mod $(PROG) 

