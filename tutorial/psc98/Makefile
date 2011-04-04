all: ;
	mpicc -I../../src/include main.c yosen.c -L../../src/lib -lestiva -lm

0: ;
	export PSC98=0; time ./a.out -solver cgs --mpi -psc98 
1: ;
	export PSC98=1; time ./a.out --mpi -psc98 
2: ;
	export PSC98=2; time ./a.out -solver cgs --mpi -psc98 
3: ;
	export PSC98=3; time ./a.out -solver cgs --mpi -psc98 
4: ;
	export PSC98=4; time ./a.out --mpi -psc98 -v -solver cgs
5: ;
	export PSC98=5; time ./a.out -solver cgs --mpi -psc98 

allclean: ;
	(cd ../../src/lib; make clean)

lib: ;
	(cd ../../src/lib/; make)