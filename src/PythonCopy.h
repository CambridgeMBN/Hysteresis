#include <cmath>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <algorithm>

#include "PythonCopy.cpp"

const double Element::mu_0 = 4 * M_PI * 1e-7;
const double Element::mu_B = 9.274e-24;
double Element::theta = M_PI / 4;
double Element::delta = 0;
std::vector<std::vector<Stack *>> Element::elementGroup; // Element list
std::vector<double>::iterator Element::exchangeIterator; // List of J_ex at interfaces
std::vector<std::pair<int, double>> Element::ddelta;
std::vector<std::vector<double>> Element::magnetisation;

double Element::H;
bool Element::isNotReversed = false;
double Element::T = 6.4; // Temperature
double Element::topLayerExchangeBias = 0; // Top layer may be exchange biased
std::vector<double> Element::exchangeList;
int Element::count = 0; // count

void Element::saveToFile() {
	std::ofstream angles;
	std::ofstream deltas;

	char fileName[100];

	double B = Element::H * Element::mu_0 * 100;
	std::cout << "B: " << B << std::endl;
	sprintf(fileName, "/Users/sriv1211/Dropbox/University of Cambridge/Term Project/hysteresis/Arrow_Data/Arrows_%i_%4.2f", (int) B, Element::T);
	angles.open(fileName);

	deltas.open("/Users/sriv1211/Dropbox/University of Cambridge/Term Project/hysteresis/delta.txt");

	std::vector<std::vector<Stack *>>::iterator list_it;
	std::vector<Stack *>::iterator stack_it;

	for (list_it = Element::elementGroup.begin();
			list_it != Element::elementGroup.end(); list_it++) {
		for (stack_it = list_it->begin(); stack_it != list_it->end();
				stack_it++) {
			angles << *(*stack_it)->phi << std::endl;

		}
	}

	std::vector<std::pair<int, double>>::iterator deltait;
	for (deltait = Element::ddelta.begin(); deltait != Element::ddelta.end(); deltait++) {
		deltas << deltait->first << "," << deltait->second << std::endl;
	}

	angles.close();
	deltas.close();

	std::cout << "Wrote to " << fileName << std::endl;
}

void Element::saveMToFile() {
	std::ofstream angles;

	char fileName[100];

	sprintf(fileName, "Arrow_Data/Mag_%4.2f", Element::T);
	angles.open(fileName);

	std::vector<std::vector<double>>::iterator list_it;
	std::vector<double>::iterator stack_it;

	for (list_it = Element::magnetisation.begin();
			list_it != Element::magnetisation.end(); list_it++) {
		angles << *list_it->begin() << ",";
		for (stack_it = list_it->begin() + 1; stack_it != list_it->end();
				stack_it++) {
			angles << *stack_it << ",";

		}
		angles << std::endl;
	}

	angles.close();
	std::cout << "Wrote to " << fileName << std::endl;
}

void Element::prepare(std::vector<double> J_ex) {
	Element::delta = 0;
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

void Element::getMagnetisation() {
	std::vector<double> magnetisation;

	magnetisation.push_back(Element::H * 0.012); // 1st element is H (in Oe)

	// Calculate individual magnetisations for the different layers
	double mTot = 0;
	for (std::vector<Stack *> stackList : Element::elementGroup) {
		double m = 0;
		for (Stack * stack : stackList) {
			m += stack->el->M(Element::T) * cos(*stack->phi);
		}
		magnetisation.push_back(m);
		mTot += m;
	}

	magnetisation.push_back(mTot); // Last element of magnetisation list is the total

	Element::magnetisation.push_back(magnetisation);
}

void Element::phaseIterate() {
	/*
	 * Iterates through all the elements in all the lists to determine the phase at a given element of a given list
	 */

	std::vector<std::vector<Stack *>>::iterator list_it;
	std::vector<Stack *>::iterator stack_it;

	list_it = Element::elementGroup.begin();
	int i = 0;

	for (list_it = Element::elementGroup.begin();
			list_it != Element::elementGroup.end(); list_it++) {
		for (stack_it = list_it->begin(); stack_it != list_it->end();
				stack_it++) {

			double H_eff = Element::H;

			Stack * prev;
			Stack * it;
			Stack * next;

			if (stack_it == list_it->begin()) {
				if (list_it == Element::elementGroup.begin()) {

					if (Element::isNotReversed) {
						H_eff += Element::topLayerExchangeBias;
					}

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

			Element::setPhase(prev, it, next, H_eff);

			i++;

			Element::count++;
		}

		Element::ddelta.push_back(std::make_pair(Element::count, Element::delta));
	}

}

void Element::createStack(int layers, double phi, Element * el) {
	std::cout << "layers: " << layers << std::endl;
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
	if (next->el->getElement().compare(it->el->getElement()) == 0) {
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


	double delta_1 = (*it->phi) - phase;

	if (fabs(delta_1) >= fabs(Element::delta)) {
//		std::cout << "count: " << Element::count << " pinit: " << (*it->phi) << " vs " << phase << " delta :" << delta_1 << std::endl;
		Element::delta = fabs(delta_1);
	}

	/*
	* the phase must change by less than pi
	*/

	if (fabs(delta_1) <= M_PI)
	{
		*it->phi = phase;
	} else {
//		std::cout << " phase: " << phase << " vs " << *it->phi << std::endl;
	}

}

void runModel() {

	new Ni(3.5);
	new Gd(7.9);
	new Ni(3.5);

	Element::setH(0.1 / Element::mu_0); // Magnetic field strength
	Element::topLayerExchangeBias = 0 / Element::mu_0;
	std::cout << "H: " << Element::getH() << std::endl;
	std::vector<double> J_ex = { -5e-22, -5e-22 };

	for (int i = 0; i <= 5000; i++) {
////	while (Element::delta >= 1e-3) {
//
		Element::prepare(J_ex);
		Element::count = 0;
		Element::phaseIterate();
//
		// Reverse array to iterate from bottom (but no top layer bias)
		Element::isNotReversed = !Element::isNotReversed;
		std::reverse(Element::elementGroup.begin(), Element::elementGroup.end());
//
		std::cout << "delta: " << Element::delta << std::endl;
	}


	Element::saveToFile();

}
