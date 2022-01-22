#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/*
double a = 5e-10;               // Lattice Const.
double b = 1e-10;               // Step Len.
double v0 = 1;                  // Potential depth eV
double E;                       // Energy eV
double dE = 0.001;              // Energy step eV
double h_bar = 6.582119569e-16; //eV.s
double m_e = 9.10938356e-31;    //e mass Kg;
double e_const = 1.602e-19;
*/
double a = 10;
double b = 8;
double v0 = 3;
double E;
double dE = 0.000001;
double h_bar = 1;
double m_e = 1;
double e_const = 1;

double cos_K_E(double, double, double, double);

int main(void)
{
    FILE *bandstrc;
    bandstrc = fopen("bandstrc.dat","w");
    for (E = -0.5; abs(E) < v0; E -= dE)
    {
        if (fabs(cos_K_E(a, b, v0, E)) <= 1)
        {
            fprintf(bandstrc, "%f\t%f\n" , E, acos(cos_K_E(a, b, v0, E))/M_PI);
        }
    }
    fclose(bandstrc);
    return 0;
}

double cos_K_E(double a, double b, double v0, double E)
{

    double alpha = sqrt(2 * m_e * fabs(E)) / h_bar;
    double beta = sqrt(2 * m_e * (v0 - fabs(E))) / h_bar;
    double d = b - a;
    double cos_ka = cos(beta * b) * cosh(alpha * d) - (beta * beta - alpha * alpha) / (2 * alpha * beta) * sin(beta * b) * sinh(alpha * d);
    //printf("%f\n",cos_ka);
    return cos_ka;
}