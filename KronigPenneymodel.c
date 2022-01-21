#include <stdio.h>
#include <math.h>

double a = 5e-10;               // Lattice Const.
double b = 1e-10;               // Step Len.
double v0 = 1;                  // Potential depth eV
double E;                       // Energy eV
double dE = 0.001;              // Energy step eV
double h_bar = 6.582119569e-16; //eV.s
double m_e = 9.10938356e-31;    //e mass Kg;
double e_const = 1.602e-19;

double K_E(double, double, double, double);
aljshalshdljhsdl
int main(void)
{
    for (E = 0; abs(E) < v0; E -= dE)
    {
        printf("%f", K_E(a, b, v0, E));
    }
    return 0;
}

double K_E(double a, double b, double v0, double E)
{

    double alpha = sqrt(2 * m_e * abs(E)) / h_bar;
    double beta = sqrt(2 * m_e * (v0 - abs(E))) / h_bar;
    double d = b - a;
    double cos_ka = cos(beta * b) * cosh(alpha * d) - (beta * beta - alpha * alpha) / (2 * alpha * beta) * sin(beta * b) * sinh(alpha * d);
    //printf("%f",cos_ka);
    if (cos_ka <= 0)
    {
        double k = acos(cos_ka) / a;
        return k;
    }
    return 0;
}