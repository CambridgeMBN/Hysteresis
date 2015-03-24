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

	Stack() {}
};

struct Magnetisation {
	double H; // external field
	std::vector<Stack> stackList;
};

class Element {
private:
	double T_c; // critical temperature
	double rho;
	double J_s;
	double static T;
	double static H;
	std::string element; // Element name, e.g. Ni, Gd
	std::vector<double> static exchangeList; // J_ex list

public:
	double static delta;
	double static const mu_0;
	double static const mu_B;
	double static theta;
	double static topLayerExchangeBias;
	bool static isNotReversed;
	int static layers;
	int static count;

	std::vector<std::vector<Stack *>> static elementGroup;
	std::vector<double>::iterator static exchangeIterator;
	std::vector<std::pair<int, double>> static ddelta;

	virtual double M(double T) = 0;
	virtual double K(double T) = 0;

	static std::vector<std::vector<double>> magnetisation;

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

	void setElement(std::string element) {
		this->element = element;
	}

	std::string getElement() {
		return element;
	}

	static void phaseIterate();

	static void setT(double T);

	static double getT();

	static void setH(double H);

	static double getH();

	static void createStack(int layers, double phi, Element * el);

	static void setPhase(Stack *prev, Stack *it, Stack *next) {
		Element::setPhase(prev, it, next, H);
	}
	static void setPhase(Stack *prev, Stack *it, Stack *next, double H_eff);

	static void getMagnetisation();

	static void prepare(std::vector<double> J_ex);

	static void saveToFile();
	static void saveMToFile();
};

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

#endif /* PYTHONCOPY_CPP2_ */
