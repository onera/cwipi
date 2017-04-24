# ------------------------------------------------------------------
#
# makefile pour fortran Intel - MacOSX
#
# Auteur : Christophe Peyret christophe@peyret.fr
# ------------------------------------------------------------------

UNAME = $(shell uname -p)

TARGET = basesNodales

### Compilers and linker ###

FC  = ifort
LD  = ifort
CC  = icc
CXX = icpc

### Paths ###

SRCDIR    = ./src
OBJDIR    = ./obj
MODS      = ./mod
SPA_EXE   = ./


### MKL ###
INCL += -I$(MKL_HOME)/include/intel64/lp64 -I${MKL_HOME}/include

#MODS +=  -module ${MKL_HOME}/include/intel64/lp64

#LIBS +=  ${MKL_HOME}/lib/libmkl_lapack95_lp64.a \
#         ${MKL_HOME}/lib/libmkl_intel_lp64.a    \
#         ${MKL_HOME}/lib/libmkl_intel_thread.a  \
#         ${MKL_HOME}/lib/libmkl_core.a -lpthread -lm

#LIBS += -lpthread -lm -mkl
LIBS += -lm -mkl

#LIBS +=  $(MKL_HOME)/lib/libmkl_lapack95_ilp64.a \
#         $(MKL_HOME)/lib/libmkl_intel_ilp64.a    \
#         $(MKL_HOME)/lib/libmkl_intel_thread.a  \
#         $(MKL_HOME)/lib/libmkl_core.a -lpthread -lm

### Compiler Flags and Link ###

FCFLAGS    = -c -O3 -fpp -mkl=sequential -arch x86_64 $(INCL) -module $(MODS)
CCFLAGS    = -c -O3                      -arch x86_64 $(INCL) -std=c99
CXXFLAGS   = -c -O3                      -arch x86_64 $(INCL) -std=c99


#Debug
#FCFLAGS   += -g -p -O0 -fpe-all=0 -traceback -ftrapuv -fp-stack-check -debug -i4 -r8 # -check all -assume realloc_lhs
#CCFLAGS   += -g
#CXXFLAGS  += -g 

LDFLAGS    =  # -cxxlib  # -openmp # -static-intel

### List of Objects ###

OBJS =  $(OBJDIR)/libmesh7.o              \
        $(OBJDIR)/table_tet_mesh.o        \
        $(OBJDIR)/modDeterminant.o        \
        $(OBJDIR)/pyramidRule.o           \
        $(OBJDIR)/baseSimplexTools.o      \
        $(OBJDIR)/baseSimplex1D.o         \
        $(OBJDIR)/baseSimplex2D.o         \
        $(OBJDIR)/baseSimplex3D.o         \
        $(OBJDIR)/basePyramid.o           \
        $(OBJDIR)/baseWedge.o             \
        $(OBJDIR)/LDLt.o                  \
        $(OBJDIR)/basesNodales.o

### Building Rules ###

$(TARGET) : $(OBJS)
	@echo "******"
	@echo Link de $(TARGET)
	@echo "*"
	$(LD) $(LDFLAGS) -o $(SPA_EXE)$(TARGET) $(OBJS) $(LIBS)
	@echo  ----------- ${TARGET} created ----------- 

compil_date :
	touch src/basesNodales.f90
	\rm -f compil_date.h
	echo "write(*,'(/a )')'------------------------------------------------------------------------'" >  compil_date.h
	echo "write(*,'( a )')'                      program Bases Nodales                             '" >> compil_date.h
	echo "write(*,'( a )')'                                                                        '" >> compil_date.h
	echo "write(*,'( a )')' Compiled on" `date` "with" `uname -n` " '"                                >> compil_date.h
	echo "write(*,'( a )')' " `svn info | grep "Changed Date"` "                                   '" >> compil_date.h
	echo "write(*,'( a )')' " `svn info | grep Revision` "                                         '" >> compil_date.h
	echo "write(*,'( a/)')'------------------------------------------------------------------------'" >> compil_date.h

clean :
	\rm -f $(OBJS) $(MODS)/*.mod

# Compilation Rules

$(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	@echo "*"
	@echo Compilation de $<
	@echo "*"
	$(CXX) $(CXXFLAGS) -o $@  $< 

$(OBJDIR)/%.o : $(SRCDIR)/%.c
	@echo "*"
	@echo Compilation de $<
	@echo "*"
	$(CC) $(CCFLAGS) -o $@  $< 

$(OBJDIR)/%.o : $(SRCDIR)/%.f
	@echo "*"
	@echo Compilation de $<
	@echo "*"
	$(FC) $(FCFLAGS) -o $@  $< 

$(OBJDIR)/%.o : $(SRCDIR)/%.f90
	@echo "*"
	@echo Compilation de $<
	@echo "*"
	$(FC) $(FCFLAGS) -o $@  $< 

$(OBJDIR)/iSend.o : $(SRCDIR)/iSend.f90 compil_date
	@echo "*"
	@echo Compilation de iSend.f90
	@echo "*"
	$(FC) $(FCFLAGS) -o $(OBJDIR)/iSend.o    $(SRCDIR)/iSend.f90

