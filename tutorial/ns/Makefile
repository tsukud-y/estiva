CC  = gcc
SRC = ../../src
SRCS = main.c A.c Rhs.c boundary_condition.c
OBJS = $(SRCS:.c=.o)
SOLVER = qmr

.c.o: 
	$(CC) -g -Wall -I${SRC}/include -c $< -o $@

all: $(OBJS)
	(cd ../../src; make)
	$(CC) $(OBJS) -g -Wall -L${SRC}/lib -L/usr/local/lib -lestiva -lm  -lefence -lmpich -lopa -lmpl -pthread

run: all
	./a.out -plotscale 0.001 -n 10 -v -solver $(SOLVER)
	(cd anime; cat *.gnuplot | gnuplot)

large: all
	./a.out -plotscale 0.001 -kn 500 -n 50 
	(cd anime; cat *.gnuplot | gnuplot)

clean: ;
	rm -f ./a.out ./a.out.core $(OBJS)

allclean: clean
	(cd ../../src; make clean)
	rm -rf ./anime
