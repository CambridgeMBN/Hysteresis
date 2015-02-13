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
#include "Energy.h"

int main() {

	double thickness = 3;
	Ni ni = Ni(thickness, M_PI, 1e-5);
	Gd gd = Gd(thickness, 0);
	double J_int = -5.e-22; // interface exchange coupling constant
	TwoLayers stack = TwoLayers(&ni, &gd, J_int);

	Energy e = Energy(stack);
//	e.iter(0.5, 10);
	e.mag(5);
	return 0;
}
