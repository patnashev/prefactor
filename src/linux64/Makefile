# Makefile for Linux

EXE       = prefactor
COMPOBJS  = precompute.o md5.o arithmetic.o edwards.o giant.o group.o lucas.o montgomery.o inputnum.o primelist.o logging.o file.o stage1.o stage2.o task.o prob.o prefactor.o
LIB_GWNUM = ../gwnum/linux64/gwnum.a

# Source directories

VPATH     = ..:../arithmetic

# insert proper dummy file to get _GWDATA data segment aligned
# on a 32 byte boundary.  Use proper dummyNN.o to do the alignment

DUMMY2 = ../gwnum/linux64/gwnum.ld
LIBS2  = -lm -lpthread -lstdc++ -static

CC  = gcc
CXX = g++
CFLAGS   = -I.. -I../gwnum                 -DX86_64 -O2 -std=c99     -Wall -Wextra
CXXFLAGS = -I.. -I../gwnum -I../arithmetic -DX86_64 -O2 -std=gnu++17 -Wall -Wextra

#########################################################################

all: $(EXE)

$(EXE): $(COMPOBJS)
	$(CC) -o $@ $(COMPOBJS) $(LIB_GWNUM) $(DUMMY2) $(LIBS2)

clean:
	rm -f $(EXE) $(COMPOBJS)
