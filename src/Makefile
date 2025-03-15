CC = g++
CFLAGS = -fPIC -c -I../include
SRC = lattice.cpp core.cpp reduction.cpp
PROGRAM = ../lib/liblat.so

target1:
	$(CC) $(CFLAGS) $(SRC)

target2:
	$(CC) -shared *.o -o $(PROGRAM)

all: target1 target2
