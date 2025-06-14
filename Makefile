all:
	cd src && make all && cd ..
	cd build && make all && cd ..
	cd test && make all && cd ..