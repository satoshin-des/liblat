CC = g++
SRC = main.cpp
LIBS = -llat
LFLAGS = -L../lib
CFLAGS = -Xlinker -rpath -Xlinker ../lib

all:
	$(CC) $(SRC) -I../include $(LFLAGS) $(LIBS) $(CFLAGS)
