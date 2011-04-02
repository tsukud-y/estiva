all:	;
	(cd /home/tsukud-y/html/estiva/tutorial/ns;make)

run: all
	./a.out -plotscale 0.001 -n 10 -v -solver cgs -adjust
	(cd anime; cat *.gnuplot | gnuplot)


allclean: ;
	(cd /home/tsukud-y/html/estiva/tutorial/ns;make allclean)