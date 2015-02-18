/*
 * PythonCopy.h
 *
 *  Created on: 15 Feb 2015
 *      Author: sriv1211
 */

#ifndef PYTHONCOPY_CPP2_
#define PYTHONCOPY_CPP2_

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

class Element;

struct Stack {
	Element * el;
	double * phi;
};

class Element {
private:
	double T_c; // critical temperature
	double rho;
	double J_s;
	double static T;
	double static H;
	int static counter;
	std::string element; // Element name, e.g. Ni, Gd
	std::vector<double> static exchangeList; // J_ex list

public:
	double static const mu_0;
	double static const mu_B;
	double static theta;
	double static topLayerExchangeBias;

	std::vector<std::vector<Stack *>> static elementGroup;
	std::vector<double>::iterator static exchangeIterator;

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

	double setElement(std::string element) {
		this->element = element;
	}

	std::string getElement() {
		return element;
	}

	static void phaseIterate();

	static void setT(double T);

	static double getT();

	static void createStack(int layers, double phi, Element * el);

	static void setPhase(Stack *prev, Stack *it, Stack *next) {
		Element::setPhase(prev, it, next, H);
	}
	static void setPhase(Stack *prev, Stack *it, Stack *next, double H_eff);

	static std::vector<double> getMagnetisation(double T);

	static void prepare(std::vector<double> J_ex);
};

// To delete
int Element::counter = 0;
const double Element::mu_0 = 4 * M_PI * 1e-7;
const double Element::mu_B = 9.274e-24;
double Element::theta = M_PI / 4;
std::vector<std::vector<Stack *>> Element::elementGroup; // Element list
//std::vector<double>::iterator Element::exchangeIterator; // List of J_ex at interfaces
//
double Element::H; // Magnetic field strength
double Element::T; // Temperature
//double Element::topLayerExchangeBias = 0; // Top layer may be exchange biased

//std::vector<double> Element::exchangeList = {};

class Gd: public Element {
private:
	double T_c = 28;
	double rho = 3.02e28;
	double J_s = 2.625e-24 * 4.45 / 2;

public:

	Gd(double thickness) :
			Element() {
		int layers = round(thickness / 0.3);
		createStack(layers, 0, this);

		setT_c(T_c);
		setRho(rho);
		setJ_s(J_s);
		setElement("Gd");
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
private:
	double T_c = 631;
	double rho = 9.14e28;
	double J_s = 1.2e-22 * 4.45 / 2;

public:

	Ni(double thickness) :
			Element() {
		int layers = round(thickness / 0.4);
		createStack(layers, M_PI, this);

		setT_c(T_c);
		setRho(rho);
		setJ_s(J_s);
		setElement("Ni");
	}

	double M(double T) {
		return 1.3 / 2; // or m?
	}

	double K(double T) {
		return 4000 / rho * 4;
	}
};

void test() {
	new Ni(7);
	new Gd(7.9);
	new Ni(7);

	std::vector<double> J_ex = {1, 2};

	Element::prepare(J_ex);
	Element::phaseIterate();

//	Element::getMagnetisation(1);
}

#endif /* PYTHONCOPY_CPP2_ */
