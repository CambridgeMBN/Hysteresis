//
// Created by Anand on 23/03/15.
//

#include "Iter.h"
//#include "PythonCopy.h"
#include <iostream>
#include <vector>
#include <tuple>


int Iterator::total = 0;
std::vector<std::vector<int>> Iterator::stacks;


int main3() {

//    std::vector<int> n = {4, 3, 2, 5};
////    int total = 0;
//    std::vector<std::vector<int>> stacks;
//    for (auto i : n) {
//
//        std::vector<int> stack;
//        for (int j = Iterator::total; j < Iterator::total + i; j++) {
//            stack.push_back(j);
//        }
//
//        stacks.push_back(stack);
//        Iterator::total += i;
//    }
//
//    Iterator it;
//    it.stacks = stacks;
//
//    for (int i = 0; i < Iterator::total; i++) {
//        int * j;
//        int * k;
//        int * l;
//        std::tie(j, k, l) = it.env();
//
//        *j += 0;
//        *k += 0;
//        *l += 0;
//
//        it.iter();
//    }
//
//    it.index = 0;
//
//    for (int i = 0; i < 3 * Iterator::total; i++) {
//
//        int * j;
//        int * k;
//        int * l;
//        std::tie(j, k, l) = it.env();
//
//        std::cout << "prev: " << *j << " it: " << *k << " next: " << *l << std::endl;
//        it.iter();
//    }

//    new Ni(3.5);
//    new Gd(7.9);
//    new Ni(3.5);
//
//    Element::setH(0.1 / Element::mu_0); // Magnetic field strength
//    Element::topLayerExchangeBias = 0 / Element::mu_0;
//    std::cout << "H: " << Element::getH() << std::endl;
//    std::vector<double> J_ex = { -5e-22, -5e-22 };
//
//    for (int i = 0; i <= 5000; i++) {
//////	while (Element::delta >= 1e-3) {
////
//        Element::prepare(J_ex);
//        Element::count = 0;
//        Element::phaseIterate();
////
//        // Reverse array to iterate from bottom (but no top layer bias)
//        Element::isNotReversed = !Element::isNotReversed;
//        std::reverse(Element::elementGroup.begin(), Element::elementGroup.end());
////
//        std::cout << "delta: " << Element::delta << std::endl;
//    }
//
//
//    Element::saveToFile();
//
//
//
//    return 0;
}