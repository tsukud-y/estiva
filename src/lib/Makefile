LIBNAME = libestiva.a
# for GNU make
# SRCS = $(shell find . -name '*.c')
# OBJS = $(subst .c,.o,$(SRCS))
SRCS = ./ary/ary.c ./esolver/esolver.c ./esolver/maxesolver.c \
	./esolver/minesolver.c ./foreach/foreach.c ./mesh/forgamma.c \
	./mesh/fp2mesh.c ./mesh/np1.c ./mesh/plt.c ./mx/clonemx.c \
	./mx/mulmx.c ./mx/mx.c ./mx/transmx.c ./op/op.c \
	./precond/precondILU.c ./precond/precondjacobi.c \
	./precond/precondnone.c ./precond/precondscaling.c ./solver/bicg.c \
	./solver/blu.c ./solver/gauss.c ./solver/pcgs.c ./solver/solver.c \
	./tmpfile/tmpfile.c ./mx/matvecmx.c
OBJS =  ./ary/ary.o ./esolver/esolver.o ./esolver/maxesolver.o \
	./esolver/minesolver.o ./foreach/foreach.o ./mesh/forgamma.o \
	./mesh/fp2mesh.o ./mesh/np1.o ./mesh/plt.o ./mx/clonemx.o \
	./mx/mulmx.o ./mx/mx.o ./mx/transmx.o ./op/op.o \
	./precond/precondILU.o ./precond/precondjacobi.o \
	./precond/precondnone.o ./precond/precondscaling.o ./solver/bicg.o \
	./solver/blu.o ./solver/gauss.o ./solver/pcgs.o ./solver/solver.o \
	./tmpfile/tmpfile.o ./mx/matvecmx.o

CFLAGS = -I../include
CPPFLAGS = -Wall -g

RM := /bin/rm
AR := /bin/ar

.SUFFIXES: .c .o
.PHONY: clean all

all: $(LIBNAME)
	

$(LIBNAME): $(OBJS)
	ar cr $(LIBNAME) $(OBJS)

.c.o: 
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

clean:
	$(RM) -f $(OBJS) $(LIBNAME)
