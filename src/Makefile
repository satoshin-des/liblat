CC = g++
CFLAGS = -fPIC -c -I../include
SRC = lattice.cpp core.cpp 

all:
	$(CC) $(CFLAGS) core.cpp -o ../build/core.o
	${CC} ${CFLAGS} lattice.cpp -o ../build/lattice.o
	${CC} ${CFLAGS} bkz.cpp -o ../build/bkz.o
	${CC} ${CFLAGS} deep_bkz.cpp -o ../build/deep_bkz.o
	${CC} ${CFLAGS} deep_lll.cpp -o ../build/deep_lll.o
	${CC} ${CFLAGS} dual_bkz.cpp -o ../build/dual_bkz.o
	${CC} ${CFLAGS} dual_deep_bkz.cpp -o ../build/dual_deep_bkz.o
	${CC} ${CFLAGS} dual_deep_lll.cpp -o ../build/dual_deep_lll.o
	${CC} ${CFLAGS} dual_enum.cpp -o ../build/dual_enum.o
	${CC} ${CFLAGS} dual_lll.cpp -o ../build/dual_lll.o
	${CC} ${CFLAGS} dual_pot_lll.cpp -o ../build/dual_pot_lll.o
	${CC} ${CFLAGS} enum.cpp -o ../build/enum.o
	${CC} ${CFLAGS} hkz.cpp -o ../build/hkz.o
	${CC} ${CFLAGS} lll.cpp -o ../build/lll.o
	${CC} ${CFLAGS} pot_bkz.cpp -o ../build/pot_bkz.o
	${CC} ${CFLAGS} pot_enum.cpp -o ../build/pot_enum.o
	${CC} ${CFLAGS} pot_lll.cpp -o ../build/pot_lll.o
	${CC} ${CFLAGS} size.cpp -o ../build/size.o

