#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>

/*
double a = 5e-10;               // Lattice Const.
double b = 1e-10;               // Step Len.
double v0 = 1;                  // Potential depth eV
double dE = 0.0001;              // Energy step eV
double h_bar = 6.582119569e-16; //eV.s
double m_e = 9.10938356e-31;    //e mass Kg;
double e_const = 1.602e-19;
*/
double a = 10;
double b = 8;
double v0 = 1.5;
double dE = 0.0001;
double h_bar = 1;
double m_e = 1;
double e_const = 1;
int bandNumber = 4;
int bandNumber2 = 4;

double cos_K_E(double, double, double, double);

int main(int argc, char *argv[])
{
    int band_i = 0;
    double ETest = 0;
    double maxEs[bandNumber];
    double minEs[bandNumber2];
    FILE *bandstrc;
    //FILE *outscf;
    bandstrc = fopen("bandstrc.dat", "w");
    //outscf = fopen("scf.dat","w");
    for (double E = 0; abs(E) <= v0; E += dE)
    {

        if (fabs(cos_K_E(a, b, v0, E)) <= 1)
        {
            if (fabs(E - ETest) > 1 * dE)
            {
                band_i += 1;
                minEs[band_i] = E;
                maxEs[band_i] = E;
            }
            ETest = E;

            if (maxEs[band_i] < E)
            {
                maxEs[band_i] = E;
            }

            if (minEs[band_i] > E)
            {
                minEs[band_i] = E;
            }

            fprintf(bandstrc, "%f\t%f\t%d\tMax:%f\tMin:%f\n", E, acos(cos_K_E(a, b, v0, E)) / M_PI, band_i,maxEs[band_i],minEs[band_i]);
    
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
    // printf("%f\n",cos_ka);
    return cos_ka;
}