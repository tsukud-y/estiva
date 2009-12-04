PROG = xmsh
SRCS = main.c
OBJS = $(SRCS:.c=.o)
CPPFLAGS = -Wall
CFLAGS = -I../../include
LDFLAGS =  -L../../lib -lestiva -lm

.PHONY: all clean

all: $(PROG)

$(PROG): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $@ $(LDFLAGS) 

.c.o:
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

clean :
	rm -f $(PROG) $(OBJS)
