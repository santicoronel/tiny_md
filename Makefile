CC       = gcc
CFLAGS	 = -O0
WFLAGS	 = -std=gnu99 -Wall -Wextra -g
LDFLAGS	 = -lm -lgomp

TARGETS		= tiny_md viz
SOURCES		= $(shell echo *.c)
OBJECTS     = core.o

all: $(TARGETS)

viz: viz.o $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) -lGL -lGLU -lglut

tiny_md: tiny_md.o $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(WFLAGS) $(CFLAGS) -c $<

clean:
	rm -f $(TARGETS) *.o *.xyz *.log .depend

.depend: $(SOURCES)
	$(CC) -MM $^ > $@

-include .depend

.PHONY: clean all
