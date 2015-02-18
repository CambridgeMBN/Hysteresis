/*
 * Energy.h
 *
 *  Created on: 12 Feb 2015
 *      Author: sriv1211
 */

#ifndef ENERGY2_H_
#define ENERGY2_H_

#include <vector>
#include <iterator>
#include <cmath>
#include <fstream>
#include <iomanip>

#include "TwoLayers.h"

class Energy {
private:

	TwoLayers stack;
	std::vector<double>::iterator phiIt;

	double theta = M_PI / 4;

	double g = 2; // g-factor
	double mu_0 = 4 * M_PI * 1e-7; // permeability
	double mu_B = 9.274e-24; // Bohr magneton

	std::vector<double *> phiTop;
	std::vector<double *> phiBottom;

public:
	Energy(TwoLayers &stack) :
		stack(stack) {
		phiTop = stack.top->getPhiList();
		phiBottom = stack.bottom->getPhiList();
	}


	void topLayer(double H, double T) {
		// start with second element
		for (auto it = phiTop.begin() + 1; it != phiTop.end(); it++) {
			std::cout << "x it: " << **it << std::endl;
		}
	}

	void verify() {
		int j = 0;
		for (auto *i: phiTop) {
			std::cout << "layer: " << j << " phi: " << *i << " \t ";
			j++;
		}

		std::cout << "---" << std::endl;

		j = 0;
		for (auto *i: phiBottom) {
			std::cout << "layer: " << j << " phi: " << i << " \t ";
			j++;
		}
	}

	void iter(double H, double T) {
//		verify();
		topLayer(H, T);

	}

//	virtual ~Energy();

};

#endif /* ENERGY2_H_ */
