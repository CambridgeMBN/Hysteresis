/*
 * Energy.h
 *
 *  Created on: 12 Feb 2015
 *      Author: sriv1211
 */

#ifndef ENERGY_H_
#define ENERGY_H_

#include <vector>
#include <iterator>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>

#include "TwoLayers.h"

class Energy {
private:

	TwoLayers stack;
	std::vector<double>::iterator phiIt;

	double theta = M_PI / 4;

	double g = 2; // g-factor
	double mu_0 = 4 * M_PI * 1e-7; // permeability
	double mu_B = 9.274e-24; // Bohr magneton
	double BFactor = g * mu_0 * mu_B; // Zeeman pre-factor

	double getIteratorEnergy(std::vector<double>::iterator it) {
	}

public:
	Energy(TwoLayers &stack) :
			stack(stack) {
	}

	double mag(double T) {
		std::vector<double> *topPhiList = stack.top->getPhiList();
		std::vector<double> *bottomPhiList = stack.bottom->getPhiList();

		/*
		 *         Mtot_Ni[k] = sum(Mt[1:nNi/2+1]*cos(phi0[1:nNi/2+1])) + \
                     sum(Mt[nNi/2+nGd+1:-0]*cos(phi0[nNi/2+nGd+1:-0]))
		 *
		 */

		double H;
		std::ofstream myfile;
		myfile.open ("hysteresis.csv");
		double x;
		int max = 40;
		double start = 4;
		H = start;
		for (int i = 0; i <= max; i++) {
			x = 4*start/max;


			double m = 0;
			iter(H, T);
			std::cout << "H: " << H << std::endl;
			for (double phase : *topPhiList) {
				m += stack.top->M(T)*cos(phase);
			}
			myfile << H << "," << m << ",";

			m = 0;
			for (double phase : *bottomPhiList) {
				m += stack.bottom->M(T)*cos(phase);
			}

			myfile << m << std::endl;


			if (i < max/2) {
							H -= x;
						} else {
							H += x;
						}

		}

		myfile.close();


	}

	void iter(double H, double T) {
		for (int i = 0; i < 1000; i++) {
			std::cout << "i: " << i << std::endl;
			double max = getPhi(H, T);
//			std::cout << "max: " << max << std::endl;
		}

		std::vector<double> *topPhiList = stack.top->getPhiList();
		std::vector<double> *bottomPhiList = stack.bottom->getPhiList();

		std::ofstream myfile;
		myfile.open ("vector.csv");
//		myfile << "Writing this to a file.\n";
		long size = topPhiList->size();
		int i = (int) size;
		i = -size;
//		int i = -1*  (int *)*topPhiList->size();
		for (double phase : *topPhiList) {
			myfile << i << "," << 5 * cos(phase) << "," << 5 * sin(phase) << std::endl;
			i++;
		}

//		i = 0;
		for (double phase : *bottomPhiList) {
			myfile << i << "," << 5 * cos(phase) << "," << 5* -sin(phase) << std::endl;
			i++;
		}

		myfile.close();

	}

	void writeToFile(double H, double T) {
		std::vector<double> *topPhiList = stack.top->getPhiList();
		std::vector<double> *bottomPhiList = stack.bottom->getPhiList();

		double B = H*Energy::mu_0;

		std::ofstream hysteresisFile;
		char fileName[64];
		sprintf(fileName, "Arrow_Data/Arrows_%4.2f%_%4.2f%", B, T);
		hysteresisFile.open (fileName);

		for (double phase : *topPhiList) {
			hysteresisFile << phase << std::endl;
		}

		for (double phase : *bottomPhiList) {
			hysteresisFile << phase << std::endl;
		}

		hysteresisFile.close();

		std::cout << " wrote to " << fileName << std::endl;
	}

	double getPhi(double H, double T) {
		Element *top = stack.top;
		Element *bottom = stack.bottom;
		double max = 0;
		double max2 = 0;

		std::vector<double> *topPhiList = stack.top->getPhiList();
		std::vector<double> *bottomPhiList = bottom->getPhiList();

		if (top->getTopBias() > 1e-6) {
			// There is a bias, start from 1st element
			phiIt = topPhiList->begin();
		} else {
			// No bias, start from 2nd element
			phiIt = topPhiList->begin() + 1;
		}

		double H_eff = H; // Effective field
		double J_int = stack.J_int;

		double JTop = top->getJ_ex(); // J for same element
		double MSatTop = top->getMSat(); // Saturation magnetisation
		double MTTop = top->M(T); // Magnetisation at finite temperature
		double KTop = top->K(T); // Anisotropic energy

		double JBottom = bottom->getJ_ex(); // J for same element
		double MSatBottom = bottom->getMSat(); // Saturation magnetisation
		double MTBottom = bottom->M(T); // Magnetisation at finite temperature
		double KBottom = bottom->K(T); // Anisotropic energy

//		std::cout << " JT: " << JBottom << " MS " << MSatBottom << " MTTop " << MTBottom << " KTop " << KBottom << std::endl;

		double ASin, ACos, ASinUp, ASinDown, ACosUp, ACosDown, ATot, BTot,
				CTot1, CTot2;
		double phase1, phase2, absPhase1, absPhase2, phase;

		int i = 0;
		std::cout << "***" << std::endl;
		for (double phase : *topPhiList) {
			std::cout << "i: " << i << " phase top: " << phase << std::endl;
			i++;
		}

		std::cout << " --- " << std::endl;
		i = 0;
		for (double phase : *bottomPhiList) {
			std::cout << "j: " << i << " phase bottom: " << phase << std::endl;
			i++;
		}
		std::cout << " *** " << std::endl;

		// Start of iterations for top stack
		for (phiIt; phiIt != topPhiList->end(); phiIt++) {
			double phiUp = *std::prev(phiIt);
			double phiDown = *std::next(phiIt);

			if (phiIt != topPhiList->begin()) {
				ASinUp = JTop * MSatTop * sin(phiUp);
				ACosUp = JTop * MSatTop * cos(phiUp);
			} else {
				H_eff = H + top->getTopBias(); // Anti-F bias in place of spin interaction with layer above
			}

			if (phiIt != topPhiList->end() - 1) {
				ASinDown = JTop * MSatTop * sin(phiDown);
				ACosDown = JTop * MSatTop * cos(phiDown);
			} else {
				phiDown = *bottomPhiList->begin();
				ASinDown = J_int * MSatBottom * sin(phiDown);
				ACosDown = J_int * MSatBottom * cos(phiDown);
			}

			ASin = ASinUp + ASinDown;
			ACos = ACosUp + ACosDown;
			ATot = sqrt(pow(ASin, 2) + pow(ACos, 2));
			BTot = BFactor * H_eff;

			CTot1 = sqrt(
					pow(MSatTop * ACos + MTTop * BTot + KTop * cos(theta), 2)
							+ pow(MSatTop * ASin + KTop * sin(theta), 2));
			CTot2 = sqrt(
					pow(MSatTop * ACos + MTTop * BTot - KTop * cos(theta), 2)
							+ pow(MSatTop * ASin - KTop * sin(theta), 2));

//			std::cout << "AS " << ASin << " AC " << ACos << " AT " << ATot << " BT " << BTot << std::endl;

			absPhase1 = asin((MSatTop * ASin + KTop * sin(theta)) / CTot1);
			absPhase2 = asin((MSatTop * ASin - KTop * sin(theta)) / CTot2);

			if (MSatTop * ACos + MTTop * BTot + KTop * cos(theta) >= 0) {
				phase1 = absPhase1;
			} else {
				phase1 = M_PI - absPhase1;
			}

			if (MSatTop * ACos + MTTop * BTot - KTop * cos(theta)) {
				phase2 = absPhase2;
			} else {
				phase2 = M_PI - absPhase2;
			}

			if (CTot1 > CTot2) {
				*phiIt = phase2;
			} else {
				*phiIt = phase1;
			}

		}

		/*
		 * BOTTOM Layer
		 */
		for (phiIt = bottomPhiList->begin(); phiIt != bottomPhiList->end() - 1;
				phiIt++) {

			double phiUp = *std::prev(phiIt);
			double phiDown = *std::next(phiIt);

			if (phiIt != bottomPhiList->begin()) {
				ASinUp = JBottom * MSatBottom * sin(phiUp);
				ACosUp = JBottom * MSatBottom * cos(phiUp);

			} else {
				phiUp = *(topPhiList->end() - 1);
				ASinUp = JBottom * MSatBottom * sin(phiUp);
				ACosUp = JBottom * MSatBottom * cos(phiUp);
//				std::cout << " exchange: " << *(topPhiList->end() -1);
			}

			if (phiIt != bottomPhiList->end() - 1) {
				ASinDown = JBottom * MSatBottom * sin(phiDown);
				ACosDown = JBottom * MSatBottom * cos(phiDown);
			}

			ASin = ASinUp + ASinDown;
			ACos = ACosUp + ACosDown;
			ATot = sqrt(pow(ASin, 2) + pow(ACos, 2));
			BTot = BFactor * H_eff;

			CTot1 = sqrt(
					pow(
							MSatBottom * ACos + MTBottom * BTot
									+ KBottom * cos(theta), 2)
							+ pow(MSatBottom * ASin + KBottom * sin(theta), 2));
			CTot2 = sqrt(
					pow(
							MSatBottom * ACos + MTBottom * BTot
									- KBottom * cos(theta), 2)
							+ pow(MSatBottom * ASin - KBottom * sin(theta), 2));

			absPhase1 = asin(
					(MSatBottom * ASin + KBottom * sin(theta)) / CTot1);
			absPhase2 = asin(
					(MSatBottom * ASin - KBottom * sin(theta)) / CTot2);

			if (MSatBottom * ACos + MTBottom * BTot + KBottom * cos(theta)
					>= 0) {
				phase1 = absPhase1;
			} else {
				phase1 = M_PI - absPhase1;
			}

			if (MSatBottom * ACos + MTBottom * BTot - KBottom * cos(theta)) {
				phase2 = absPhase2;
			} else {
				phase2 = M_PI - absPhase2;
			}

			double init = *phiIt;
			if (CTot1 > CTot2) {
				*phiIt = phase2;
			} else {
				*phiIt = phase1;
			}
		}

		i = 0;
		for (double phase : *topPhiList) {
			std::cout << "i: " << i << " phase top: " << phase << std::endl;
			i++;
		}

		std::cout << " --- " << std::endl;
		i = 0;
		for (double phase : *bottomPhiList) {
			std::cout << "j: " << i << " phase bottom: " << phase << std::endl;
			i++;
		}

		std::cout << "max : " << max << " max2: " << max2 << std::endl;

		writeToFile(H, T);
	}

};

#endif /* ENERGY_H_ */
