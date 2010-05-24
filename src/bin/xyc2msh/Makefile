PROG =  xyc2msh
SRCS = Delaunay.c  FILE.c  oldestiva.c  stack.c  xyc2msh.c
OBJS = $(SRCS:.c=.o)

CPPFLAGS = -Wall -g -O2
CFLAGS = -I../../include -fno-strict-aliasing
LDFLAGS = -L../../lib -lestiva -lm

.PHONY: all clean

all: $(PROG)

$(PROG) : $(OBJS)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OBJS) -o $@ $(LDFLAGS)

.c.o:
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

clean :
	rm -f $(OBJS) $(PROG)
