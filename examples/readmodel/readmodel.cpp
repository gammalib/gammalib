#include <iostream>
#include "GammaLib.hpp"
int main(void) {
    GModels models("source.xml");
    std::cout << models << std::endl;
    return 0;
}

