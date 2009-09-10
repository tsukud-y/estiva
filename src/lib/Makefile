all: ;
	cc -I../include -c op/op.c
	cc -I../include -c ary/ary.c
	cc -I../include -c mx/mx.c
	cc -I../include -c solver/solver.c
	cc -I../include -c solver/gauss.c
	cc -I../include -c solver/blu.c
	cc -I../include -c solver/pcgs.c
	cc -I../include -c esolver/esolver.c
	cc -I../include -c esolver/maxesolver.c
	cc -I../include -c esolver/minesolver.c
	ar cr libestiva.a *.o


clean: ;
	rm -f libestiva.a *.o
