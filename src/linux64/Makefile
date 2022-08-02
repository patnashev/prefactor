# Makefile for Linux

EXE       = prefactor
EXE_F     = prefactor_f
LIB_GWNUM = ../gwnum/linux64/gwnum.a

COMPOBJS_COMMON = md5.o arithmetic.o edwards.o giant.o lucas.o montgomery.o poly.o inputnum.o primelist.o logging.o file.o stage1.o stage2.o stage2poly.o task.o prob.o params.o
COMPOBJS   = $(COMPOBJS_COMMON) group.o   prefactor.o
COMPOBJS_F = $(COMPOBJS_COMMON) group_f.o prefactor_f.o factoring.o

# Source directories

VPATH     = ..:../arithmetic

# insert proper dummy file to get _GWDATA data segment aligned
# on a 32 byte boundary.  Use proper dummyNN.o to do the alignment

DUMMY2 = ../gwnum/linux64/gwnum.ld
LIBS2  = -lm -lpthread -lstdc++ -static

CC  = gcc
CXX = g++
CFLAGS   = -I.. -I../gwnum                           -DX86_64 -O2 -std=c99     -Wall -Wextra
CXXFLAGS = -I.. -I../gwnum -I../arithmetic -DEXT_DAC -DX86_64 -O2 -std=gnu++17 -Wall -Wextra

#########################################################################

all: $(EXE)

factoring: $(EXE_F)

$(EXE): $(COMPOBJS)
	$(CC) -o $@ $(COMPOBJS) $(LIB_GWNUM) $(DUMMY2) $(LIBS2)

$(EXE_F): $(COMPOBJS_F)
	$(CC) -o $@ $(COMPOBJS_F) $(LIB_GWNUM) $(DUMMY2) $(LIBS2)

prefactor_f.o: prefactor.cpp
	$(CXX) -o $@ -c $(CXXFLAGS) -DFACTORING $<

group_f.o: group.cpp
	$(CXX) -o $@ -c $(CXXFLAGS) -UEXT_DAC $<

clean:
	rm -f $(EXE) $(COMPOBJS)
