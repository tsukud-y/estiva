CC  = gcc
CFLAGS = -g -Wall
SRCDIR = ../../src
INCLUDEDIR = -I$(SRCDIR)/include
LIBDIR = -L$(SRCDIR)/lib -L/usr/local/lib -L/opt/mpich2/lib
SRCS = main.c A.c Rhs.c boundary_condition.c Re.c cylindermesh.c
OBJS = $(SRCS:.c=.o)
LIBS =  -lestiva -lm -lefence -lmpich -lopa -lmpl -pthread

.c.o: 
	$(CC) -c $(CFLAGS) $(INCLUDEDIR) $< -o $@

all: $(OBJS)
	(cd ../../src; make)
	$(CC)    $(CFLAGS) $(OBJS) $(LIBDIR) $(LIBS)

clean: ;
	rm -f ./a.out ./a.out.core $(OBJS)

allclean: clean
	(cd ../../src; make clean)
	rm -rf ./anime
