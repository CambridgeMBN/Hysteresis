/*
 * PythonCopy.h
 *
 *  Created on: 15 Feb 2015
 *      Author: sriv1211
 */

#ifndef PYTHONCOPY_H_
#define PYTHONCOPY_H_

#include <iostream>
#include <cmath>
#include <vector>
#include <iterator>
#include <iomanip>      // std::setprecision

class Element;

struct Stack;

class Element {
private:
	double T_c; // critical temperature
	double rho;
	double J_s;

public:

	double static const mu_0;
	double static const mu_B;
	double static theta;
	std::vector<double> phiList;

	virtual double M(double T) = 0;
	virtual double K(double T) = 0;

	Element() {

	}

	void setT_c(double T_c) {
		this->T_c = T_c;
	}

	double getT_c() {
		return T_c;
	}

	void setRho(double rho) {
		this->rho = rho;
	}

	double getRho() {
		return rho;
	}

	void setJ_s(double J_s) {
		this->J_s = J_s;
	}

	double getJ_s() {
		return J_s;
	}

	void createStack(int layers, double phi) {
		for (int i = 0; i < layers; i++) {
			phiList.push_back(phi);
		}
	}
};

const double Element::mu_0 = 4 * M_PI * 1e-7;
const double Element::mu_B = 9.274e-24;
double Element::theta = M_PI / 4;

class Gd: public Element {
public:
	double T_c = 28;
	double rho = 3.02e28;
	double J_s = 2.625e-22 * 4.45 / 2;

	Gd(double thickness, double phi) :
			Element() {
		int layers = round(thickness / 0.3);
		createStack(layers, phi);

		setT_c(T_c);
		setRho(rho);
		setJ_s(J_s);
	}

	double M(double T) {
		double MGd_min = 0.3414;
		double m = (7 - MGd_min) / 2 * (1 + tanh(0.1143 * (T_c - T))) + MGd_min;
		return m / 2; // or m?
	}

	double K(double T) {
		double k = 17500 / rho * pow((M(T) / 7), 3) * 4;
		return k;
	}
};

class Ni: public Element {
public:
	double T_c = 631;
	double rho = 9.14e28;
	double J_s = 1.2e-22 * 4.45 / 2;

	Ni(double thickness, double phi) :
			Element() {
		int layers = round(thickness / 0.4);
		createStack(layers, phi);

		setT_c(T_c);
		setRho(rho);
		setJ_s(J_s);
	}

	double M(double T) {
		return 1.3 / 2; // or m?
	}

	double K(double T) {
		return 4000 / rho * 4;
	}
};

struct Stack {

	Element * top;
	Element * down;
	double J_ex;

	void magneticEnergy(double B, double T) {
		double H = B / Element::mu_0;
		double theta = Element::theta;

		double MS_up, J_up, phiUp;

		std::cout << " down: " << down->M(0) << std::endl;

		for (int j = 0; j < 1; j++) {
			std::cout << " " << std::endl;
			std::vector<double>::iterator it = top->phiList.begin();

			for (int i = 1; i < 19; i++) {

				if (it == top->phiList.begin()) {
					MS_up = 0;
					J_up = 0;
					phiUp = 0;
				} else {
					phiUp = *std::prev(it);
					MS_up = top->M(0);
					J_up = top->getJ_s();
				}

				double K = top->K(T);
				double M_0 = top->M(0);
				double M_T = top->M(T);

				double MS_down, J_down, phiDown;
				if (it == top->phiList.end() - 1) {
					std::cout << "end " << std::endl;
					MS_down = down->M(0);
					J_down = J_ex;
					phiDown = *down->phiList.begin();
				} else {
					MS_down = top->M(0);
					J_down = top->getJ_s();
					phiDown = *std::next(it);
				}

				/*
				 *         print 'i ', i, ' H ', H, ' Ku ', Ku, ' theta ', theta, ' Mt ', Mt, ' Ms ', Ms, ' MSup', Msup,
				 *         ' phiup ', phiup,
				 *         ' Jup ', Jup, ' MsDn ', Msdn, ' PhiDn ', phidn, ' JDn ', Jdn
				 */

				std::cout << " i " << i << " H " << H << " Ku " << K << " theta " <<  theta << " MT " << M_T << " MS " << M_0 <<
						" phiUP " << phiUp << " Jup " << J_up << " MSDn " << MS_down <<
							" phiDN " << phiDown << " JDn " << J_down << std::endl;

				double ASin = J_up * MS_up * sin(phiUp)
						+ J_down * MS_down * sin(phiDown);
				double ACos = J_up * MS_up * cos(phiUp)
						+ J_down * MS_down * cos(phiDown);
				double BTot = Element::mu_0 * 2 * Element::mu_B * H;

				double CTot1 = sqrt(
						pow(M_0 * ACos + M_T * BTot + K * cos(theta), 2)
								+ pow(M_0 * ASin + K * sin(theta), 2));
				double CTot2 = sqrt(
						pow(M_0 * ACos + M_T * BTot - K * cos(theta), 2)
								+ pow(M_0 * ASin - K * sin(theta), 2));

				double absPhase1 = asin((M_0 * ASin + K * sin(theta)) / CTot1);
				double absPhase2 = asin((M_0 * ASin - K * sin(theta)) / CTot2);

				double phase1 =
						(M_0 * ACos + M_T * BTot + K * cos(theta) >= 0) ?
								absPhase1 : M_PI - absPhase1;
				double phase2 =
						(M_0 * ACos + M_T * BTot - K * cos(theta) >= 0) ?
								absPhase2 : M_PI - absPhase2;

				double phase = (CTot1 < CTot2) ? phase2 : phase1;

				//	        print 'a ', Asin, ' c ', Acos, ' ', Atot, ' Bt ', Btot, ' Ctot1 ', Ctot1, ' abs1 ', absphase1, ' abs2 ',
				// absphase2, ' ct1 ', Ctot2, \
				 ' p1 ', phase1, ' p2 ', phase2, ' p ', phase
				//             print 'i ', i, 'k ', Ku, ' abs1 ', absphase1, ' abs2 ', absphase2,
				// ' ct1 ', Ctot2, ' ct2 ', Ctot2, ' p1 ', phase1, ' p2 ', phase2, ' p ', phase
				if (i > 17) {
					std::cout << std::setprecision(15) << "i: " << i << " k " << K <<
						" abs1 " << absPhase1 << " abs2 " << absPhase2 << " ct1 " << CTot1 << " ct2 " << CTot2 <<
						" p1 " << phase1 << " p2 " << phase2 << " p " << phase
							<< std::endl;

				}

				*it = phase;
				it++;
			}
		}

	}
};

void test() {
	Ni * ni = new Ni(7, M_PI / 2);
	Gd * gd = new Gd(7.9, 0);

	Stack s;
	s.top = ni;
	s.down = gd;
	s.J_ex = -5e-22;

	s.magneticEnergy(-0.5, 6.4);
}

#endif /* PYTHONCOPY_H_ */
