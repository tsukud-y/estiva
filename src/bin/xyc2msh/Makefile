PROG =  xyc2msh
SRCS = Delaunay.c  FILE.c  oldestiva.c  stack.c  xyc2msh.c
OBJS = Delaunay.o  FILE.o  oldestiva.o  stack.o  xyc2msh.o

CPPFLAGS = -Wall
CFLAGS = -I../../include
LDFLAGS = -L../../lib -lestiva -lm

.PHONY: all clean

all: $(PROG)

$(PROG) : $(OBJS)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OBJS) -o $@ $(LDFLAGS)

.c.o:
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

clean :
	rm -f $(OBJS) $(PROG)
