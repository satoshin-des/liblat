CC = g++
CFLAGS = -fPIC -c
SRC = lattice.cpp
PROGRAM = ../lib/liblat.so

target1:
	$(CC) $(CFLAGS) $(SRC) -I../include -o ../lib/lattice.o

target2:
	$(CC) -shared ../lib/*.o -o $(PROGRAM)

all: target1 target2
