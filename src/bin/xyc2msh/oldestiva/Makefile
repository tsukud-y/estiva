default : all

all : estiva.o confary.o matprop.o solver.o numgrid.o op.o matutil.o

estiva.o : estiva.h estiva.c
	gcc -c estiva.c

confary.o : estiva.h confary.h confary.c
	gcc -c confary.c

matprop.o : matprop.c
	gcc -c matprop.c

solver.o :  solver.c 
	gcc -c solver.c

numgrid.o : numgrid.h numgrid.c
	gcc -c numgrid.c

op.o : op.c
	gcc -c op.c

matutil.o : matutil.c
	gcc -c matutil.c	

clean : ;
	rm -f *.o 
