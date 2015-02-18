#include <cmath>
#include <iostream>
#include <fstream>
#include <stdio.h>

#include "PythonCopy.cpp"

const double Element::mu_0 = 4 * M_PI * 1e-7;
const double Element::mu_B = 9.274e-24;
double Element::theta = M_PI / 4;
std::vector<std::vector<Stack *>> Element::elementGroup; // Element list
std::vector<double>::iterator Element::exchangeIterator; // List of J_ex at interfaces

double Element::H = 0.1 / Element::mu_0; // Magnetic field strength
double Element::T = 6.4; // Temperature
double Element::topLayerExchangeBias = 0; // Top layer may be exchange biased
std::vector<double> Element::exchangeList;

void Element::saveToFile() {
	std::ofstream angles;

	char fileName [100];

	double B = Element::H * Element::mu_0 * 100;
	std::cout << "B: " << B << std::endl;
	sprintf(fileName, "Arrow_Data/Arrows_%i_%4.2f",
			(int) B, Element::T);
	angles.open(fileName);

	std::vector<std::vector<Stack *>>::iterator list_it;
	std::vector<Stack *>::iterator stack_it;

	for (list_it = Element::elementGroup.begin();
			list_it != Element::elementGroup.end(); list_it++) {
		for (stack_it = list_it->begin(); stack_it != list_it->end();
				stack_it++) {
			angles << *(*stack_it)->phi << std::endl;

		}

	}

	angles.close();
}

void Element::prepare(std::vector<double> J_ex) {
	Element::exchangeList = J_ex;
	Element::exchangeIterator = Element::exchangeList.begin();
}

void Element::setH(double H) {
	Element::H = H;
}

double Element::getH() {
	return Element::H;
}

void Element::setT(double T) {
	Element::T = T;
}

double Element::getT() {
	return Element::T;
}

std::vector<double> Element::getMagnetisation(double T) {
	std::vector<double> magnetisation;

	// Calculate individual magnetisations for the different layers
	for (std::vector<Stack *> stackList : Element::elementGroup) {
		double m = 0;
		for (Stack * stack : stackList) {
			m += stack->el->M(T) * cos(*stack->phi);
		}
		magnetisation.push_back(m);
	}

	// Calculate total magnetisation
	double total = 0;
	for (double m : magnetisation) {
		total += m;
	}

	magnetisation.push_back(total); // Last element of magnetisation list is the total

	return magnetisation;
}

void Element::phaseIterate() {
	/*
	 * Iterates through all the elements in all the lists to determine the phase at a given element of a given list
	 */

	std::vector<std::vector<Stack *>>::iterator list_it;
	std::vector<Stack *>::iterator stack_it;

	for (list_it = Element::elementGroup.begin();
			list_it != Element::elementGroup.end(); list_it++) {
		for (stack_it = list_it->begin(); stack_it != list_it->end();
				stack_it++) {

			Stack * prev;
			Stack * it;
			Stack * next;

			if (stack_it == list_it->begin()) {
				if (list_it == Element::elementGroup.begin()) {
					prev = *stack_it; // first element of first list is mirrored
				} else {
					prev = *(std::prev(list_it)->end() - 1); // last element of previous list
				}
			} else {
				prev = *std::prev(stack_it); // previous element of same stack_it list
			}

			it = *stack_it;

			if (stack_it == list_it->end() - 1) {
				if (list_it == Element::elementGroup.end() - 1) {
					next = it; // last element is mirrored as next
				} else {
					next = (*std::next(list_it)->begin()); // first element of next list
				}
			} else {
				next = *std::next(stack_it); // next element of same stack_it list
			}

			Element::setPhase(prev, it, next);

//			std::cout << " prev: " << *prev->phi << " of "
//					<< prev->el->getElement() << " \t it: " << *it->phi
//					<< " of " << it->el->getElement() << " \t next "
//					<< *next->phi << " of " << next->el->getElement()
//					<< std::endl;
		}
	}
}

void Element::createStack(int layers, double phi, Element * el) {
	std::vector<Stack *> stackList;

	for (int i = 0; i < layers; i++) {
		Stack *s = new Stack();
		s->el = el;
		s->phi = new double(phi);
		stackList.push_back(s);
	}

	Element::elementGroup.push_back(stackList);
}

void Element::setPhase(Stack *prev, Stack *it, Stack *next, double H_eff) {
	double K = it->el->K(Element::T);
	double M_0 = it->el->M(0);
	double M_T = it->el->M(Element::T);

	double MS_up = prev->el->M(0);
	double phiUp = *prev->phi;
	double J_up;

	/*
	 * J coupling top upper layer depends if they it is the same element (done by comparing string for element name). If not,
	 * get value from exchangeIterator pointer, then iterate it.
	 */
	if (prev->el->getElement().compare(it->el->getElement()) == 0) {
		J_up = prev->el->getJ_s();
	} else {
		J_up = *Element::exchangeIterator;
	}

	double MS_down = next->el->M(0);
	double phiDown = *next->phi;
	double J_down;

	// Same for bottom layer
	if (prev->el->getElement().compare(it->el->getElement()) == 0) {
		J_down = next->el->getJ_s();
	} else {
		J_down = *Element::exchangeIterator;
		Element::exchangeIterator++;
	}

	double ASin = J_up * MS_up * sin(phiUp) + J_down * MS_down * sin(phiDown);
	double ACos = J_up * MS_up * cos(phiUp) + J_down * MS_down * cos(phiDown);
	double BTot = Element::mu_0 * 2 * Element::mu_B * H_eff;

	double CTot1 = sqrt(
			pow(M_0 * ACos + M_T * BTot + K * cos(Element::theta), 2)
					+ pow(M_0 * ASin + K * sin(Element::theta), 2));
	double CTot2 = sqrt(
			pow(M_0 * ACos + M_T * BTot - K * cos(Element::theta), 2)
					+ pow(M_0 * ASin - K * sin(Element::theta), 2));

	double absPhase1 = asin((M_0 * ASin + K * sin(Element::theta)) / CTot1);
	double absPhase2 = asin((M_0 * ASin - K * sin(Element::theta)) / CTot2);

	double phase1 =
			(M_0 * ACos + M_T * BTot + K * cos(Element::theta) >= 0) ?
					absPhase1 : M_PI - absPhase1;
	double phase2 =
			(M_0 * ACos + M_T * BTot - K * cos(Element::theta) >= 0) ?
					absPhase2 : M_PI - absPhase2;

	double phase = (CTot1 < CTot2) ? phase2 : phase1;

	*it->phi = phase;
}

void test() {
	new Ni(7);
	new Gd(7.9);
//	new Ni(7);

	std::cout << "H: " << Element::getH() << std::endl;
	std::vector<double> J_ex = { 1, 2 };



	for (int i = 0; i < 100; i++) {
		Element::prepare(J_ex);
		Element::phaseIterate();
	}

	Element::saveToFile();

//	Element::getMagnetisation(1);
}
