# Unix

CC   = cc
F77  = f77
#F77  = g77 -fno-second-underscore
OBJ  = o
EXE  =
FOUT = -o fmain
RM   = rm

# Windows

#CC   = cl
#F77  = f77
#OBJ  = obj
#EXE  = .exe
#FOUT =
#RM   = del

# Fortran to C interface

#F2C = -DUPPERCASE
#F2C = -DUPPERCASE_
#F2C = -DUPPERCASE__
#F2C = -DLOWERCASE
#F2C = -DLOWERCASE_
#F2C = -DLOWERCASE__

# additional Fortran libraries

FLIBS =

# targets

test : fmain$(EXE)
	fmain$(EXE)

fmain$(EXE) : fmain.f csub.$(OBJ)
	$(F77) $(FOUT) fmain.f csub.$(OBJ) $(FLIBS)

csub.$(OBJ) : csub.c Makefile
	$(CC) $(F2C) -c csub.c

clean :
	-$(RM) *.$(OBJ)
	-$(RM) fmain$(EXE)
