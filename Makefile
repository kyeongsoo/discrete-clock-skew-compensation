# --- Compiler and Flags ---
CC = gcc
CFLAGS = -Wall -O3
# CFLAGS = -Wall -g
LDFLAGS = -lm
LDLIBS = -lmpfr -lgmp

SRCS = csc.c csc-adds-ana.c csc-adds-ana2.c csc-ana.c drand48.c
OBJS = $(SRCS:.c=.o)
SHARED_OBJS = csc.o

# macro definitions
CFLAGS += -DCSC_INT_SIZE=4
# CFLAGS += -DCSC_INT_SIZE=8
# CFLAGS += -DCSC_INT_SIZE=16
CFLAGS += -DDEBUG=1

all: csc-adds-ana csc-ana csc-analysis csc-debug1 csc-debug2 csc-debug3

csc-adds-ana: $(SHARED_OBJS) drand48.o csc-adds-ana.o
	$(CC) -o $@ $^ $(LDFLAGS) $(LDLIBS)

csc-ana: $(SHARED_OBJS) drand48.o csc-ana.o
	$(CC) -o $@ $^ $(LDFLAGS) $(LDLIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# A phony target for cleaning up generated files
.PHONY: all clean

clean:
	rm -f csc csc-adds-ana csc-adds-ana2 csc-ana $(OBJS)
