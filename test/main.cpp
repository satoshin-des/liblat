#include <iostream>

#include "core.h"
#include "lattice.h"

int main()
{
    Lattice<int> lat(50, 50);
    lat.setRandom(40, 40, 1000, 10000);
    // lat.setRandom(40, 40, 1000, 10000);
    std::cout << lat;
    lat.computeGSO();
    std::vector<int> v = lat.mulVecBasis(lat.ENUM(1000));
    print(v);
    // lat.dualBKZ(20, 0.99);
    lat.deepBKZ(30, 0.99);
    // lat.dualDeepBKZ(30, 0.99);
    // lat.BKZ(20, 0.99);
    // lat.BKZ(6, 0.99);
    std::cout << lat.b1Norm() << std::endl;
    // std::cout << lat.volume(false) << std::endl;

    return 0;
}
