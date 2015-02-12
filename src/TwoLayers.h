/*
 * TwoLayers.h
 *
 *  Created on: 12 Feb 2015
 *      Author: sriv1211
 */

#ifndef TWOLAYERS_H_
#define TWOLAYERS_H_

#include "Element.h"

struct TwoLayers {
	Element *top;
	Element *bottom;
	double J_int; // interface J

	TwoLayers(Element *top, Element *bottom, double J_int) : top(top), bottom(bottom), J_int(J_int)
	{}
};


#endif /* TWOLAYERS_H_ */
