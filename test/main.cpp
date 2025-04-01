#include <iostream>

#include "core.h"
#include "lattice.h"

int main()
{
    Lattice<int> lat(40, 40);
    lat.setGoldesteinMayerLattice(10000, 1000);
    // lat.setRandom(40, 40, 1000, 10000);
    std::cout << lat;

    std::vector<int> v = lat.enumShortVec();
    print(v);
    lat.LLL(0.99);
    // lat.deepBKZ(30, 0.99);
    // lat.BKZ(20, 0.99);
    // lat.BKZ(30, 0.99);
    std::cout << lat;

    return 0;
}
