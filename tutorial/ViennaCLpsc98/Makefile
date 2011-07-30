a.out: main.cpp
	cc -c yosen.c
	c++ main.cpp yosen.o -lOpenCL -lm

run: a.out
	export PSC98=1; time ./a.out
	export PSC98=2; time ./a.out
	export PSC98=3; time ./a.out
	export PSC98=4; time ./a.out
	export PSC98=5; time ./a.out

clean: ;
	rm -f a.out *~  *.o
