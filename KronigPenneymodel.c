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
double a = 8.90;
double b = 8.0;
double v0 = 1.5;
double dE = 0.005;
double h_bar = 1;
double m_e = 1;
double e_const = 1;
int bandNumber = 4;

double cos_K_E(double a, double b, double v0, double E)
{
    double alpha = sqrt(2 * m_e * fabs(E)) / h_bar;
    double beta = sqrt(2 * m_e * (v0 - fabs(E))) / h_bar;
    double d = b - a;
    double cos_ka = cos(beta * b) * cosh(alpha * d) - (beta * beta - alpha * alpha) / (2 * alpha * beta) * sin(beta * b) * sinh(alpha * d);
    // printf("%f\n",cos_ka);
    return cos_ka;
}
void printMatrix(int m, int n, double matrix[m][n])
{
    int i, j;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            printf("%lf\t", matrix[i][j]);
        }
        printf("\n");
    }
}

int main(int argc, char *argv[])
{
    int k_grid = v0 / dE;
    int k_points[bandNumber];
    double E_k[bandNumber][k_grid][2];
    int band_i = 0;
    double ETest = 0;
    double maxEs[bandNumber];
    double minEs[bandNumber];

    for (double E = 0; E <= v0; E += dE)
    {

        if (fabs(cos_K_E(a, b, v0, E)) <= 1)
        {
            if (fabs(E - ETest) > 2 * dE)
            {   
                band_i++;
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
            E_k[band_i][k_points[band_i]][0] = E;
            E_k[band_i][k_points[band_i]][1] = acos(cos_K_E(a, b, v0, E)) / M_PI;
            k_points[band_i]++;
        }
    }

    // print to file bandstrc.dat
    FILE *bandstrc;
    bandstrc = fopen("bandstrc.dat", "w");
    fprintf(bandstrc, "E\tk\tBand#\n");
    for (int i = 1; i <= band_i; i++)
    {
        for (int k = 0; k < k_points[i]; k += 1)
        {
            fprintf(bandstrc, "%f\t%f\t%d\n", E_k[i][k][0], E_k[i][k][1], i);
        }
    }
    fclose(bandstrc);
    // print to file bandstrc.dat

    FILE *gap;
    gap = fopen("gap.dat", "w");
    fprintf(gap, "Band #\tE_max\tE_min\n");
    for (int i = 1; i <= band_i; i++)
    {
        fprintf(gap, "%d\t%f\t%f\n", i, maxEs[i], minEs[i]);
    }
    
    return 0;
}