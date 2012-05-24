SRCS    = src/main.cpp
PROG    = bin/a.out
LDFLAGS = -Iinclude -Llib -lestivaplus -lestiva -lOpenCL -lf2c -lm
CFLAGS  = -Wall -o $(PROG)

$(PROG): $(SRCS)
	(cd lib; make)
	c++ $(SRCS) $(LDFLAGS) $(CFLAGS)

run: $(PROG)
	$(PROG)

clean: ;
	rm -f $(PROG)
	(cd lib; make clean)
