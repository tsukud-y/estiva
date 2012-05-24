all: ;
	cc  -c -g -Wall *.c   -I../include
	c++ -c -g -Wall *.cpp -I../include
	ar cr libestivaplus.a *.o
	rm *.o

clean: ;
	rm -f libestivaplus.a