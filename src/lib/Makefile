all: ;
	cc -I../include -c -g  -Wall op/op.c
	cc -I../include -c -g  -Wall ary/ary.c
	cc -I../include -c -g  -Wall mx/mx.c
	cc -I../include -c -g  -Wall mx/transmx.c
	cc -I../include -c -g  -Wall mx/mulmx.c
	cc -I../include -c -g  -Wall precond/precondnone.c
	cc -I../include -c -g  -Wall precond/precondscaling.c
	cc -I../include -c -g  -Wall solver/solver.c
	cc -I../include -c -g  -Wall solver/gauss.c
	cc -I../include -c -g  -Wall solver/blu.c
	cc -I../include -c -g  -Wall solver/pcgs.c
	cc -I../include -c -g  -Wall solver/bicg.c
	cc -I../include -c -g  -Wall esolver/esolver.c
	cc -I../include -c -g  -Wall esolver/maxesolver.c
	cc -I../include -c -g  -Wall esolver/minesolver.c
	cc -I../include -c -g  -Wall mesh/fp2mesh.c
	cc -I../include -c -g  -Wall mesh/plt.c
	cc -I../include -c -g  -Wall mesh/forgamma.c
	cc -I../include -c -g  -Wall mesh/np1.c
	cc -I../include -c -g  -Wall foreach/foreach.c
	cc -I../include -c -g  -Wall tmpfile/tmpfile.c
	ar cr libestiva.a *.o


clean: ;
	rm -f libestiva.a *.o
