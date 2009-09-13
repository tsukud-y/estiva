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
	cc -I../include -c mesh/fp2mesh.c
	cc -I../include -c mesh/plt.c
	cc -I../include -c mesh/pltuv.c
	cc -I../include -c mesh/forgamma.c
	cc -I../include -c mesh/np1.c
	cc -I../include -c foreach/foreach.c
	cc -I../include -c tmpfile/tmpfile.c
	ar cr libestiva.a *.o


clean: ;
	rm -f libestiva.a *.o
