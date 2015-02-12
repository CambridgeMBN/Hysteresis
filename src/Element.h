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
	double topBias; // Top layer may be AF biased (e.g. by FeMn)
	std::vector<double> phiList;

	/*
	 * Check this over
	 */
//	double thicknessWeight; // Not sure how this relates to thickness

public:
	virtual double M(double T) = 0;
	virtual double K(double T) = 0;

	void makeElementStack(int layers, double phi)
	{
		for (int i = 0; i < layers; i++) {
			phiList.push_back(phi);
		}

		std::cout << "Layers: " << layers << " phi: " << phi << std::endl;
	}

	double getJ_ex() {
		return J_ex;
	}

	Element(double topBias) : topBias(topBias) {}

};

class Ni : public Element {
private:
	double rho = 9.14e28;
	double J_ex = 1.2e-22*4.45;
	double thicknessWeight = 0.4;

public:

	Ni(double thickness, double phi) : Ni(thickness, phi, 0) {}

	Ni(double thickness, double phi, double topBias) : Element(topBias)
	{
		int layers = round(thickness/thicknessWeight);
		this->makeElementStack(layers, phi);
	}

	double M(double T)
	{
		return 1.3; // Assumed constant as T_c is large
	}

	double K(double T)
	{
		return 4e3/rho*2;
	}

};

class Gd : public Element {
private:
	double rho = 3.02e28;
	double J_ex = 2.625e-24*4.45;

	double T_c = 28;
	double M_min = 0.3414;

	double thicknessWeight = 0.3;

public:

	Gd(double thickness, double phi) : Gd(thickness, phi, 0) {}

	Gd(double thickness, double phi, double topBias) : Element(topBias)
	{
		int layers = round(thickness/thicknessWeight);
		this->makeElementStack(layers, phi);
	}

	void setLayers(double thickness) {
			std::cout << "weight: " << this->thicknessWeight << std::endl;
		}

	double M(double T) {
		// (7.-MGd_min)/2*(1+tanh(0.1143*(Tc[1]-T)))+MGd_min
		double m = (7 - M_min)/2*(1 + tanh(0.1143*(T_c - T)))+M_min;
		return m;
	}

	double K(double T) {
		// 17500./N_Gd*(MT[1]/7*2)**3/4])*3/1.5
		double k = 17500/rho * pow( (M(T)/14), 3/4) * 2;
		return k;
	}

};
#endif /* ELEMENT_H_ */
