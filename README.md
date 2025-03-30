# liblat

This is a C++ library for lattice reduction.

## How to Install

First, you need to clone this repository to your computer:

```shell
$ git clone https://github.com/satoshin-des/liblat.git
```

Use ``cd`` command to change the directories to liblat:

```shell
$ cd liblat
```

Next, use ``make`` command to compile the source codes:

```shell
$ make
cd src && make all && cd ..
make[1]: Entering directory 'hoge/liblat/src'
g++ -fPIC -c -I../include lattice.cpp core.cpp reduction.cpp
g++ -shared *.o -o ../lib/liblat.so
make[1]: Leaving directory 'hoge/liblat/src'
cd test && make all && cd ..
make[1]: Entering directory 'hoget/liblat/test'
g++ main.cpp -L../lib -llat -I../include -Xlinker -rpath -Xlinker ../lib
make[1]: Leaving directory 'hoge/liblat/test'
```
