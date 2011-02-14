SRC = ../../src

all: ;
	cc *.c -I${SRC}/include -L${SRC}/lib -lestiva -lm
	./a.out -kn 60 -velocity 100.0 -plotscale 0.005	
	(cd anime; cat *.gnuplot | gnuplot)

clean: ;
	rm -f ./a.out *.o
