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
	./tmpfile/tmpfile.c ./mx/matvecmx.c ./mx/psolvemx.c ./solver/cgs.c \
	./solver/bicgstab.c ./fgetline/fgetline.c ./fgetline/chomp.c \
	./que/que.c ./fgetline/fsize.c ./que/forq.c

OBJS = $(SRCS:.c=.o)
CFLAGS = -I../include
CPPFLAGS = -Wall -g -O2

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
