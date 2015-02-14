/*
 * Element.h
 *
 *  Created on: 12 Feb 2015
 *      Author: sriv1211
 */

#ifndef ELEMENT_H_
#define ELEMENT_H_

#import <math.h>
#import <string>
#import <vector>

class Element {
private:
	double J_ex; // exchange energy of same spin
	double rho; // density
	double topBias = 0; // Top layer may be AF biased (e.g. by FeMn)
	double MSat; // Saturation magnetisation
	std::vector<double> phiList;

	/*
	 * Check this over
	 */
//	double thicknessWeight; // Not sure how this relates to thickness
public:
	virtual double M(double T) = 0;
	virtual double K(double T) = 0;
	virtual double getMSat() = 0;
	virtual double getJ_ex() = 0;
	static double count;

	void makeElementStack(int layers, double phi) {
		for (int i = 0; i < layers; i++) {
			phiList.push_back(phi); // should be phi
		}

		count += layers;
		std::cout << "Layers: " << layers << " phi: " << phi << std::endl;
	}

	double getTopBias() {
		return topBias;
	}

	std::vector<double> * getPhiList() {
		return &phiList;
	}

	Element(double topBias) :
			topBias(topBias) {
	}

};

double Element::count = 0;

class Ni: public Element {
private:
	double rho = 9.14e28;
	double J_ex = 1.2e-22 * 4.45;
	double thicknessWeight = 0.4;
	double MSat = 1.3;

public:

	Ni(double thickness, double phi) :
			Ni(thickness, phi, 0) {
	}

	Ni(double thickness, double phi, double topBias) :
			Element(topBias) {
		int layers = round(thickness / thicknessWeight);
		this->makeElementStack(layers, phi);

		std::cout << "Ni K: " << this->K(5) << " M(T) " << this->M(5) << " , " << 4e3/rho*2*2 << std::endl;
	}

	double M(double T) {
		return 1.3; // Assumed constant as T_c is large
	}

	double K(double T) {
		return 4e3 / rho * 2;
	}

	double getMSat() {
		return MSat;
	}

	double getJ_ex() {
		return J_ex;
	}

};

class Gd: public Element {
private:
	double rho = 3.02e28;
	double J_ex = 2.625e-24 * 4.45;

	double T_c = 28;
	double M_min = 0.3414;
	double MSat = 7;

	double thicknessWeight = 0.3;

public:

	Gd(double thickness, double phi) :
			Gd(thickness, phi, 0) {
	}

	Gd(double thickness, double phi, double topBias) :
			Element(topBias) {
		int layers = round(thickness / thicknessWeight);
		this->makeElementStack(layers, phi);

		std::cout << "Gd K: " << this->K(5) << " M(T) " << this->M(5) << std::endl;
	}

	void setLayers(double thickness) {
		std::cout << "weight: " << this->thicknessWeight << std::endl;
	}

	double M(double T) {
		// (7.-MGd_min)/2*(1+tanh(0.1143*(Tc[1]-T)))+MGd_min
		double m = (7 - M_min) / 2 * (1 + tanh(0.1143 * (T_c - T))) + M_min;
		return m;
	}

	double K(double T) {
		// 17500./N_Gd*(MT[1]/7*2)**3/4])*3/1.5
		double k = 17500 / rho * pow((M(T) / 14), 3 / 4) * 2;
		return k;
	}

	double getMSat() {
		return MSat;
	}

	double getJ_ex() {
		return J_ex;
	}

};
#endif /* ELEMENT_H_ */

///*
// * Element.h
// *
// *  Created on: 12 Feb 2015
// *      Author: sriv1211
// */
//
//#ifndef ELEMENT_H_
//#define ELEMENT_H_
//
//#import <math.h>
//#import <string>
//#import <vector>
//
//class Element {
//private:
//	double J_ex; // exchange energy of same spin
//	double rho; // density
//	double topBias = 0; // Top layer may be AF biased (e.g. by FeMn)
//	double MSat; // Saturation magnetisation
//	std::vector<double *> phiList;
//
//	/*
//	 * Check this over
//	 */
////	double thicknessWeight; // Not sure how this relates to thickness
//public:
//	virtual double M(double T) = 0;
//	virtual double K(double T) = 0;
//	virtual double getMSat() = 0;
//	virtual double getJ_ex() = 0;
//
//	void makeTopElementStack(int layers, double phi) {
//
//		for (auto it = phiList.begin(); it != phiList.end(); it++) {
//			std::cout << "1st it: " << **it << std::endl;
//		}
//
//		double * top = new double(phi);
//		phiList.push_back(top);
//		phiList.push_back(top);
//		for (int i = 2; i < layers; i++) {
//			phiList.push_back(new double(phi)); // should be phi
//		}
//
//		for (auto it = phiList.begin(); it != phiList.end(); it++) {
//			std::cout << "2nd it: " << **it << std::endl;
//		}
//
//	}
//
//	void makeBottomElementStack(int layers, double phi) {
//		for (int i = 0; i < layers -2; i++) {
//			phiList.push_back(new double(phi)); // should be phi
//		}
//
//		double * bottom = new double(phi);
//		phiList.push_back(bottom);
//		phiList.push_back(bottom);
//
//	}
//
//
//
//	double getTopBias() {
//		return topBias;
//	}
//
//	std::vector<double *>  getPhiList() {
//		return phiList;
//	}
//
//	Element(double topBias) :
//			topBias(topBias) {
//	}
//
//};
//
//class Ni: public Element {
//private:
//	double rho = 9.14e28;
//	double J_ex = 1.2e-22 * 4.45;
//	double thicknessWeight = 0.4;
//	double MSat = 1.3;
//
//public:
//
//	Ni(double thickness, double phi) :
//			Ni(thickness, phi, 0) {
//	}
//
//	Ni(double thickness, double phi, double topBias) :
//			Element(topBias) {
//		int layers = round(thickness / thicknessWeight);
//		this->makeTopElementStack(layers, phi);
//	}
//
//	double M(double T) {
//		return 1.3; // Assumed constant as T_c is large
//	}
//
//	double K(double T) {
//		return 4e3 / rho * 2;
//	}
//
//	double getMSat() {
//		return MSat;
//	}
//
//	double getJ_ex() {
//		return J_ex;
//	}
//
//};
//
//class Gd: public Element {
//private:
//	double rho = 3.02e28;
//	double J_ex = 2.625e-24 * 4.45;
//
//	double T_c = 28;
//	double M_min = 0.3414;
//	double MSat = 7;
//
//	double thicknessWeight = 0.3;
//
//public:
//
//	Gd(double thickness, double phi) :
//			Gd(thickness, phi, 0) {
//	}
//
//	Gd(double thickness, double phi, double topBias) :
//			Element(topBias) {
//		int layers = round(thickness / thicknessWeight);
//		this->makeBottomElementStack(layers, phi);
//	}
//
//	double M(double T) {
//		// (7.-MGd_min)/2*(1+tanh(0.1143*(Tc[1]-T)))+MGd_min
//		double m = (7 - M_min) / 2 * (1 + tanh(0.1143 * (T_c - T))) + M_min;
//		return m;
//	}
//
//	double K(double T) {
//		// 17500./N_Gd*(MT[1]/7*2)**3/4])*3/1.5
//		double k = 17500 / rho * pow((M(T) / 14), 3 / 4) * 2;
//		return k;
//	}
//
//	double getMSat() {
//		return MSat;
//	}
//
//	double getJ_ex() {
//		return J_ex;
//	}
//
//};
//#endif /* ELEMENT_H_ */
