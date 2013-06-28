include ../Makefile.in

BINDIR = ../Graphs

INCLUDES = -I../Lib 
CFLAGS = $(COPTIONS) $(OPTFLAGS) $(INCLUDES)

LIBS = -lmetis -lm
LD = $(CC) $(LDOPTIONS) -L. -L..
METISLIB = ../libmetis.a



MTESTOBJS = mtest.o ../Programs/io.o 


.c.o:
	$(CC) $(CFLAGS) -c $*.c

default: $(BINDIR)/mtest

$(BINDIR)/mtest: $(MTESTOBJS) $(METISLIB)
	$(LD)  -o $@ $(MTESTOBJS) $(LIBS)
	chmod 755 $@

clean:
	rm -f *.o

realclean:
	rm -f *.o ;\
        rm -f $(BINDIR)/mtest 
