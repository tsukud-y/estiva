SRC = ../../src

SRCS = main.c A.c Rhs.c boundary_condition.c
OBJS = $(SRCS:.c=.o)

.c.o: 
	$(CC) -Wall -I${SRC}/include -c $< -o $@

all: $(OBJS)
	(cd ../../src; make)
	$(CC) $(OBJS) -Wall -L${SRC}/lib -lestiva -lm

run: all
	./a.out -plotscale 0.005
	(cd anime; cat *.gnuplot | gnuplot)

clean: ;
	rm -f ./a.out ./a.out.core $(OBJS)

allclean: clean
	(cd ../../src; make clean)
	rm -rf ./anime