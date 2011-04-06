SOLVER = cg
V = 

yosen: ;
	(cd ../../src/lib/; make)
	mpicc -I../../src/include main.c yosen.c  -L../../src/lib -lestiva -lm 

honsen: ;
	(cd ../../src/lib/; make)
	mpicc -I../../src/include main.c honsen.c  -L../../src/lib -lestiva -lm

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
	(cd ../../src/lib; make clean)
	rm -f ./a.out ./a.out.core


