#include <iostream>

#include "core.h"
#include "lattice.h"

int main()
{
    Lattice<int> lat(40, 40);
    lat.setRandom(40, 40, 1000, 10000);
    lat.computeGSO();
    std::cout << lat;

    std::vector<int> v = lat.enumShortVec();
    print(v);
    // lat.LLL(0.99);
    lat.HKZ(0.99);
    // lat.BKZ(30, 0.99);
    std::cout << lat;

    return 0;
}
