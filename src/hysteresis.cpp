//============================================================================
// Name        : hysteresis.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <cmath>
#include <vector>
#include <iterator>

#include "Element.h"
#include "TwoLayers.h"

int main() {

	double thickness = 2.5;
	Gd gd = Gd(thickness, 0);
	Ni ni = Ni(thickness, M_PI);
	double J_int = -5.e-22; // interface exchange couplign constant
	TwoLayers stack = TwoLayers(&gd, &ni, J_int);



//	LayerList top = LayerList(&gd);
//	LayerList bottom = LayerList(&ni);
//	TwoLayers stack = TwoLayers(top, bottom);


//	std::vector<int> top, bottom;
//	std::vector<int>::iterator it;
//	bool isBiased = true;
//
//	for (int i = 0; i < 5; i++) {
//		top.push_back(i);
//	}
//
//	for (int i = 5; i < 10; i++) {
//		bottom.push_back(i);
//	}
//
//	it = (isBiased) ? top.begin() : top.begin() + 1; // set start based on bias
//
//	for (it; it != top.end(); it++) {
//
//		auto upLayer = std::prev(it);
//		auto downLayer = std::next(it);
//
//		if (it != top.begin()) {
//			std::cout << "prev: " << *upLayer;
//		} else {
//			std::cout << "bias";
//		}
//
//		std::cout << " i: " << *it;
//
//		if (it != top.end() - 1) {
//			std::cout << " next: " << *downLayer;
//		} else {
//			std::cout << " exchange: " << *bottom.begin();
//		}
//
//		 std::cout << std::endl;
//
//	}
//
//	std::cout << " --- " << std::endl;
//
//	for (it = bottom.begin(); it != bottom.end() -1; it++) {
//
//		auto prev = std::prev(it);
//		auto next = std::next(it);
//
//		if (it != bottom.begin()) {
//			std::cout << "prev: " << *prev;
//		} else {
//			std::cout << " exchange: " << *(top.end() -1);
//		}
//
//		std::cout << " i: " << *it;
//
//		if (it != bottom.end() - 1) {
//			std::cout << " next: " << *next << std::endl;
//		}
//
//	}

	return 0;
}
