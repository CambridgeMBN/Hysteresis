//
// Created by Anand on 23/03/15.
//

#include "Iter.h"

#include <iostream>
#include <vector>
#include <tuple>

int main() {

    std::vector<int> stack = {5, 3, 2};

    int total = 0;
    for (auto i : stack) {
        total += i;
    }

    std::cout << "t: " << total << std::endl;

    int index = 3;

    bool indexFound = false;
    while (!indexFound) {

        

        indexFound = true;
    }


    return 0;
}