#include "Data.h"
#include "Vector.h"

#ifndef min
template <class T> static inline T min(T x, T y) {return (x<y)?x:y;}
#endif

#ifndef max 
template <class T> static inline T max(T x, T y) {return (x>y)?x:y;}
#endif

double CVector::WDotXt(const double *w, const CFeature_node *x_t)
{
    double f_t = 0.0;
    while(x_t->index != -1)
    {
        f_t += w[x_t->index-1]*x_t->value;
        x_t++;
    }
    return f_t;
}

void CVector::WAddXt(double *w, const CFeature_node *x_t,double y_t)
{    
    while(x_t->index != -1)
    {
        w[x_t->index-1] += y_t*x_t->value;
        x_t++;
    }    
}


//x_t .* x_t * sigma 
double CVector::V3_dot(const double *s,const CFeature_node *x_t)
{
    double v_t = 0.0;
    while(x_t->index != -1)
    {
        v_t += x_t->value*x_t->value*s[x_t->index-1];
        x_t++;
    }
    return v_t;
}

double CVector::Norm2(const CFeature_node *x_t)
{
    double norm = 0.0;
    while(x_t->index != -1)
    {
        norm += (x_t->value)*(x_t->value);
        x_t++;
    }
    norm = sqrt(norm);
    return norm;
}
