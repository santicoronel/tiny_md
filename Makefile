CC      = gcc
CFLAGS	= -g -march=native -O3 -ffast-math -fopenmp
WFLAGS	= -std=c11 -Wall -Wextra -Werror
LDFLAGS	= -lm

TARGETS	= tiny_md viz
SOURCES	= $(shell echo *.c)
OBJECTS = core.o wtime.o

tiny_md: tiny_md.o $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

all: $(TARGETS)

viz: viz.o $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) -lGL -lGLU -lglut


%.o: %.c
	$(CC) $(WFLAGS) $(CPPFLAGS) $(CFLAGS) -c $<

clean:
	rm -f $(TARGETS) *.o *.xyz *.log .depend

.depend: $(SOURCES)
	$(CC) -MM $^ > $@

-include .depend

.PHONY: clean all