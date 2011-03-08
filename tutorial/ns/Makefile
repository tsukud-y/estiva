SRC = ../../src

all: ;
	cc *.c -Wall -I${SRC}/include -L${SRC}/lib -lestiva -lm

run: all
	./a.out -kn 100 -velocity 100.0 -plotscale 0.005	
	(cd anime; cat *.gnuplot | gnuplot)

clean: ;
	rm -f ./a.out *.o
