#include <iostream>

#include "core.h"
#include "lattice.h"

int main()
{
    Lattice<int> lat(10, 10); // 10-demensional full-rank lattice
    // set as a random lattice
    lat.setRandom(10, 10, 1000, 10000);

    // print the lattice basis
    std::cout << lat;

    // compute GSO-information
    lat.computeGSO();

    // compute the shortest vector
    std::vector<int> v = lat.mulVecBasis(lat.ENUM(1000));
    print(v);

    // lattice basis reduction
    lat.LLL(0.99);
    // lat.dualBKZ(20, 0.99);
    //lat.deepBKZ(30, 0.99);
    // lat.dualDeepBKZ(30, 0.99);
    // lat.BKZ(20, 0.99);
    // lat.BKZ(6, 0.99);
    std::cout << lat.b1Norm() << std::endl;
    // std::cout << lat.volume(false) << std::endl;

    return 0;
}
