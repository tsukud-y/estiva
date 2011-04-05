SOLVER = qmr

all: ;
	mpicc -I../../src/include main.c yosen.c -L../../src/lib -lestiva -lm 

0: ;
	export PSC98=0; time ./a.out -solver $(SOLVER) --mpi -psc98 
1: ;
	export PSC98=1; time ./a.out -solver $(SOLVER) --mpi -psc98
2: ;
	export PSC98=2; time ./a.out -solver $(SOLVER) --mpi -psc98 
3: ;
	export PSC98=3; time ./a.out -solver $(SOLVER) --mpi -psc98 
4: ;
	export PSC98=4; time ./a.out -solver $(SOLVER) --mpi -psc98
5: ;
	export PSC98=5; time ./a.out -solver $(SOLVER) --mpi -psc98

run: 1 2 3 4 5


allclean: ;
	(cd ../../src/lib; make clean)

lib: ;
	(cd ../../src/lib/; make)
