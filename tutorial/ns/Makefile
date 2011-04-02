SRC = ../../src

SRCS = main.c A.c Rhs.c boundary_condition.c
OBJS = $(SRCS:.c=.o)

.c.o: 
	$(CC) -Wall -I${SRC}/include -c $< -o $@

all: $(OBJS)
	(cd ../../src; make)
	$(CC) $(OBJS) -Wall -L${SRC}/lib -L/usr/local/lib -lestiva -lm -lefence

run: all
	./a.out -plotscale 0.001 -n 10 -v -solver cgs
	(cd anime; cat *.gnuplot | gnuplot)

large: all
	./a.out -plotscale 0.001 -kn 500 -n 50 
	(cd anime; cat *.gnuplot | gnuplot)

clean: ;
	rm -f ./a.out ./a.out.core $(OBJS)

allclean: clean
	(cd ../../src; make clean)
	rm -rf ./anime