PROG = unestiva
# for GNU make
# SRCS = $(shell find . -name '*.c')
# OBJS = $(subst .c,.o,$(SRCS))
SRCS = main.c
OBJS = $(SRCS:.c=.o)


.PHONY: all clean

all: $(PROG)

$(PROG): $(OBJS)
	$(CC) $(OBJS) -o $@

.c.o: 
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $<

clean: ;
	rm -f $(OBJS) $(PROG)
