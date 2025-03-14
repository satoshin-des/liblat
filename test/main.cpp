#include <iostream>

#include "lattice.h"

int main()
{
    Lattice<int> lat;
    lat.setRandom(10, 10, 1000, 10000);
    lat.computeGSO();
    std::cout << lat;

    return 0;
}
