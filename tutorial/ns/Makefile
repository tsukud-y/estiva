all: ;
	cc *.c -I${HOME}/include -L${HOME}/lib -lestiva -lm
	./a.out -kn 60 -velocity 100.0 -plotscale 0.005	
	(cd anime; cc *.c; ./a.out)