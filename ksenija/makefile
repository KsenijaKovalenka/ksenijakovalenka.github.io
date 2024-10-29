.SUFFIXES: .f90
################# Compiler flag #####################
FC = mpiifort     # Fortran compiler

################ Link flags (and libraries) #########
PARPACKDIR = ../ARPACK_PARPACK-master/ARPACK_MINE/ARPACK_MINE/Lib
#MKLROOT = /opt/intel/oneapi/mkl/2024.0


FCFLAGS = -c  -O2 -xcore-avx512 -ipo -qopt-zmm-usage=high  #from makefile2
FLFLAGS =
INC = -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include
LIBS = $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a \
       $(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a -L$(MKLROOT)/lib/intel64 \
                -lmkl_intel_lp64 -lmkl_intel_thread \
                -lmkl_core -liomp5 -lpthread -lm -ldl -lstdc++ \
         $(PARPACKDIR)/myparpack.a

################ Work directory ######################
WKDIR = .

################ Program files ######################
MODUL = $(WKDIR)/DEBUG_MODULE.f90 $(WKDIR)/PCONTEXT_MODULE.f90 $(WKDIR)/STAT_MODULE.f90
MAIN = $(WKDIR)/main.f90 $(WKDIR)/minimisation.f90 $(WKDIR)/utilities.f90

# Source files and object files
SRC = $(MODUL) $(MAIN)
OBJ = $(SRC:.f90=.o) 

################ Execution #########################
# Rule to compile .f90 into .o files
%.o: %.f90
	$(FC) $(FCFLAGS) -o $@ $<


############# Create Exe file #######################
PROGRAM = myparpack
$(PROGRAM): $(OBJ)
	$(FC) $(INC) $(FLFLAGS) -o $@ $^ $(LIBS) 
	chmod ugo+r *
	chmod ugo+x $(PROGRAM)

############# Make tar file #########################
tar:
	rm -f FSURFER.tar
	tar -cvf FSURFER.tar $(SRC) $(PROGRAM)

############ Clean ##################################
clean:
	rm -f $(WKDIR)/*.mod $(WKDIR)/*.o $(PROGRAM)

