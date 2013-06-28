include ../Makefile.in

BINDIR = ..

INCLUDES = -I../Lib 
CFLAGS = $(COPTIONS) $(OPTFLAGS) $(INCLUDES)

LIBS = -lmetis -lm
LD = $(CC) $(LDOPTIONS) -L. -L..
METISLIB = ../libmetis.a



PMETISOBJS = pmetis.o io.o 
KMETISOBJS = kmetis.o io.o 
OEMETISOBJS = oemetis.o io.o smbfactor.o
ONMETISOBJS = onmetis.o io.o smbfactor.o
MESH2DUALOBJ = mesh2dual.o io.o 
MESH2NODALOBJ = mesh2nodal.o io.o 
PARTDMESHOBJ = partdmesh.o io.o 
PARTNMESHOBJ = partnmesh.o io.o 
GRAPHCHKOBJ = graphchk.o io.o 


.c.o:
	$(CC) $(CFLAGS) -c $*.c

default: $(BINDIR)/pmetis$(VERNUM) $(BINDIR)/kmetis$(VERNUM) \
         $(BINDIR)/oemetis$(VERNUM) $(BINDIR)/onmetis$(VERNUM) $(BINDIR)/mesh2dual$(VERNUM) \
         $(BINDIR)/mesh2nodal$(VERNUM) $(BINDIR)/partdmesh$(VERNUM) $(BINDIR)/partnmesh$(VERNUM) \
         $(BINDIR)/graphchk$(VERNUM)

$(BINDIR)/pmetis$(VERNUM): $(PMETISOBJS) $(METISLIB)
	$(LD)  -o $@ $(PMETISOBJS) $(LIBS)
	chmod 755 $@

$(BINDIR)/kmetis$(VERNUM): $(KMETISOBJS) $(METISLIB)
	$(LD)  -o $@ $(KMETISOBJS) $(LIBS)
	chmod 755 $@

$(BINDIR)/oemetis$(VERNUM): $(OEMETISOBJS) $(METISLIB)
	$(LD)  -o $@ $(OEMETISOBJS) $(LIBS)
	chmod 755 $@

$(BINDIR)/onmetis$(VERNUM): $(ONMETISOBJS) $(METISLIB)
	$(LD)  -o $@ $(ONMETISOBJS) $(LIBS)
	chmod 755 $@

$(BINDIR)/mesh2dual$(VERNUM): $(MESH2DUALOBJ) $(METISLIB)
	$(LD)  -o $@ $(MESH2DUALOBJ) $(LIBS)
	chmod 755 $@

$(BINDIR)/mesh2nodal$(VERNUM): $(MESH2NODALOBJ) $(METISLIB)
	$(LD)  -o $@ $(MESH2NODALOBJ) $(LIBS)
	chmod 755 $@

$(BINDIR)/partdmesh$(VERNUM): $(PARTDMESHOBJ) $(METISLIB)
	$(LD)  -o $@ $(PARTDMESHOBJ) $(LIBS)
	chmod 755 $@

$(BINDIR)/partnmesh$(VERNUM): $(PARTNMESHOBJ) $(METISLIB)
	$(LD)  -o $@ $(PARTNMESHOBJ) $(LIBS)
	chmod 755 $@

$(BINDIR)/graphchk$(VERNUM): $(GRAPHCHKOBJ) $(METISLIB)
	$(LD)  -o $@ $(GRAPHCHKOBJ) $(LIBS)
	chmod 755 $@


clean:
	rm -f *.o

realclean:
	rm -f *.o ;\
        rm -f $(BINDIR)/pmetis$(VERNUM) ;\
        rm -f $(BINDIR)/kmetis$(VERNUM) ;\
        rm -f $(BINDIR)/oemetis$(VERNUM) ;\
        rm -f $(BINDIR)/onmetis$(VERNUM) ;\
        rm -f $(BINDIR)/mesh2dual$(VERNUM) ;\
        rm -f $(BINDIR)/mesh2nodal$(VERNUM) ;\
        rm -f $(BINDIR)/partdmesh$(VERNUM) ;\
        rm -f $(BINDIR)/partnmesh$(VERNUM) ;\
        rm -f $(BINDIR)/graphchk$(VERNUM)
