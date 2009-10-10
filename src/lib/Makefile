all: ;
	cc -I../include -c -Wall op/op.c
	cc -I../include -c -Wall ary/ary.c
	cc -I../include -c -Wall mx/mx.c
	cc -I../include -c -Wall solver/solver.c
	cc -I../include -c -Wall solver/gauss.c
	cc -I../include -c -Wall solver/blu.c
	cc -I../include -c -Wall solver/pcgs.c
	cc -I../include -c -Wall solver/bicg.c
	cc -I../include -c -Wall esolver/esolver.c
	cc -I../include -c -Wall esolver/maxesolver.c
	cc -I../include -c -Wall esolver/minesolver.c
	cc -I../include -c -Wall mesh/fp2mesh.c
	cc -I../include -c -Wall mesh/plt.c
	cc -I../include -c -Wall mesh/forgamma.c
	cc -I../include -c -Wall mesh/np1.c
	cc -I../include -c -Wall foreach/foreach.c
	cc -I../include -c -Wall tmpfile/tmpfile.c
	ar cr libestiva.a *.o


clean: ;
	rm -f libestiva.a *.o
