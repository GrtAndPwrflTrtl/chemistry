OB_INC=../openbabel-2.3.2/include/
OB_LIB_DIR=.
OB_LIB=openbabel

IDIR =./
CC=g++
OPT= -g -pg #-O1
CFLAGS= $(OPT) -I$(IDIR) -I$(OB_INC) -l$(OB_LIB) -lpthread
#
ODIR=./obj

DEPS = EdgeAggregator.h \
	EdgeAnnotation.h \
	HyperEdge.h \
	HyperEdgeMultiMap.h \
	HyperGraph.h \
	HyperNode.h \
	Instantiator.h \
	Linker.h \
	Molecule.h \
	PebblerHyperEdge.h \
	PebblerHyperGraph.h \
	PebblerHyperNode.h \
	Rigid.h \
	Utilities.h \
	Atom.h \
	Bond.h \
	AtomT.h \
	IdFactory.h \
	OBWriter.h \
	Options.h \
	Validator.h \
	obgen.h \
	Constants.h \
	Thread_Pool.h

_OBJ = Atom.o \
	Instantiator.o \
	Linker.o\
	Main.o \
	Molecule.o \
	Rigid.o \
	Utilities.o \
	Atom.o \
	Bond.o \
	AtomT.o \
	OBWriter.o \
	IdFactory.o \
	Options.o \
	Validator.o \
	obgen.o \
	Constants.o

OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) 

synth: $(OBJ)
	$(CC) $^ $(CFLAGS) -o $@

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core synth.exe synth.exe.stackdump $(INCDIR)/*~
