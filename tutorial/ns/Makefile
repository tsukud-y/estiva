SRC = ../../src

SRCS = main.c A.c Rhs.c boundary_condition.c
OBJS = $(SRCS:.c=.o)

.c.o: 
	$(CC) -Wall -I${SRC}/include -c $< -o $@

all: $(OBJS)
	(cd ../../src; make)
	$(CC) $(OBJS) -Wall -L${SRC}/lib -lestiva -lm

run: all
	./a.out -kn 100 -velocity 100.0 -plotscale 0.005	
	(cd anime; cat *.gnuplot | gnuplot)

clean: ;
	rm -f ./a.out $(OBJS)
