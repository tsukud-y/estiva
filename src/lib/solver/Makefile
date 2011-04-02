all:	;
	(cd ../../..//tutorial/ns;make)

run: all
	./a.out -plotscale 0.001 -n 10 -v -solver cgs -adjust
	(cd anime; cat *.gnuplot | gnuplot)


allclean: ;
	(cd ../../..//tutorial/ns;make allclean)