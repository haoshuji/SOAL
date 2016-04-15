#ifndef VECTOR_H
#define VECTOR_H

#include "Data.h"

class CVector{
public:
	static double WDotXt(const double *w, const CFeature_node *x_t);
	static void WAddXt(double *w, const CFeature_node *x_t, double y_t);
	static double V3_dot(const double *s, const CFeature_node *x_t);
	static double Norm2(const CFeature_node *x_t);
};

#endif
