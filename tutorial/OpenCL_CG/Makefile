# -lduma
SOLVER = cg
V = 

yosen: ../../src/lib/libestiva.a
	gcc -g  -I../../src/include main.c yosen.c tocuda.c  cg.c Rf4.c Rf5.c -L../../src/lib -lestiva -lm -lOpenCL

honsen: ../../src/lib/libestiva.a
	gcc -g -I../../src/include main.c honsen.c tocuda.c cg.c Rf4.c  Rf5.c -L../../src/lib -lestiva -lm -lOpenCL

../../src/lib/libestiva.a: ;
	(cd ../../src; make)

clean: a.out
	rm -f a.out *~

0: ;
	export PSC98=0; time ./a.out -solver $(SOLVER) $(V)
1: ;
	export PSC98=1; time ./a.out -solver $(SOLVER) $(V)
2: ;
	export PSC98=2; time ./a.out -solver $(SOLVER) $(V)
3: ;
	export PSC98=3; time ./a.out -solver $(SOLVER) $(V)
4: ;
	export PSC98=4; time ./a.out -solver $(SOLVER) $(V)
5: ;
	export PSC98=5; time ./a.out -solver $(SOLVER) $(V)

run: 1 2 3 4 5


allclean: ;
	rm -f ./a.out ./a.out.core


