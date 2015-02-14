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
#include "Energy.h"
#include "TwoLayers.h"

int main() {

//	std::vector<double *> x;
//	double	y = 5;
//	x.push_back(&y);
//	x.push_back(&y);
//
//	for (double i = 0.0; i < 3; i++) {
//		double * r = new double(i);
//		x.push_back(r);
//	}
//
//	for (auto i : x) {
//		std::cout << " i(x) " << i << std::endl;
//	}
//
//	std::cout << "y: " << &y << " x " << x[0];


	double thickness = 2.5;
	Ni ni = Ni(thickness, M_PI, 1e-5);
	Gd gd = Gd(thickness, 0);
	double J_int = -5.e-22; // interface exchange coupling constant
	TwoLayers stack = TwoLayers(&ni, &gd, J_int);

	Energy e = Energy(stack);
	e.iter(0.5, 10);
////	e.mag(5);
	return 0;
}
