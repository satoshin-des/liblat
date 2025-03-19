#include <iostream>

#include "core.h"
#include "lattice.h"

int main()
{
    Lattice<int> lat(10, 10);
    lat.setRandom(10, 10, 1000, 10000);
    lat.computeGSO();
    std::cout << lat;

    std::vector<long> v = lat.enumShortVec();
    print(v);

    lat.deepLLL(0.99);
    std::cout << lat;

    return 0;
}
