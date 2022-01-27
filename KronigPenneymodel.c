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
int bandNumber = 5;

double cos_K_E(double a, double b, double v0, double E)
{
    double alpha = sqrt(2 * m_e * fabs(E)) / h_bar;
    double beta = sqrt(2 * m_e * (v0 - fabs(E))) / h_bar;
    double d = b - a;
    double cos_ka = cos(beta * b) * cosh(alpha * d) - (beta * beta - alpha * alpha) / (2 * alpha * beta) * sin(beta * b) * sinh(alpha * d);
    // printf("%f\n",cos_ka);
    return cos_ka;
}
void gaussEliminationLS(int m, int n, double a[m][n], double x[n - 1])
{
    int i, j, k;
    for (i = 0; i < m - 1; i++)
    {
        // Partial Pivoting
        for (k = i + 1; k < m; k++)
        {
            // If diagonal element(absolute vallue) is smaller than any of the terms below it
            if (fabs(a[i][i]) < fabs(a[k][i]))
            {
                // Swap the rows
                for (j = 0; j < n; j++)
                {
                    double temp;
                    temp = a[i][j];
                    a[i][j] = a[k][j];
                    a[k][j] = temp;
                }
            }
        }
        // Begin Gauss Elimination
        for (k = i + 1; k < m; k++)
        {
            double term = a[k][i] / a[i][i];
            for (j = 0; j < n; j++)
            {
                a[k][j] = a[k][j] - term * a[i][j];
            }
        }
    }
    // Begin Back-substitution
    for (i = m - 1; i >= 0; i--)
    {
        x[i] = a[i][n - 1];
        for (j = i + 1; j < n - 1; j++)
        {
            x[i] = x[i] - a[i][j] * x[j];
        }
        x[i] = x[i] / a[i][i];
    }
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
void fit(int k_point, double x[k_point], double y[k_point])
{
    int i, j;

    // the degree of polynomial to be used:
    int n = 3;
    // an array of size 2*n+1 for storing N, Sig xi, Sig xi^2, ...., etc. which are the independent components of the normal matrix
    double X[2 * n + 1];
    for (i = 0; i <= 2 * n; i++)
    {
        X[i] = 0;
        for (j = 0; j < k_point; j++)
        {
            X[i] = X[i] + pow(x[j], i);
        }
    }
    // the normal augmented matrix
    double B[n + 1][n + 2];
    // rhs
    double Y[n + 1];
    for (i = 0; i <= n; i++)
    {
        Y[i] = 0;
        for (j = 0; j < k_point; j++)
        {
            Y[i] = Y[i] + pow(x[j], i) * y[j];
        }
    }
    for (i = 0; i <= n; i++)
    {
        for (j = 0; j <= n; j++)
        {
            B[i][j] = X[i + j];
        }
    }
    for (i = 0; i <= n; i++)
    {
        B[i][n + 1] = Y[i];
    }
    double A[n + 1];
    printf("The polynomial fit is given by the equation:\n");
    //printMatrix(n + 1, n + 2, B);
    gaussEliminationLS(n + 1, n + 2, B, A);
    for (i = 0; i <= n; i++)
    {
        printf("%lfx^%d+", A[i], i);
    }
    printf("...\n");
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
    for (int i = 1; i <= band_i; i++)
    {
        printf("fiting");

        double x[k_points[i]];
        double y[k_points[i]];
        for (int k = 0; k < k_points[i]; k++)
        {
            x[k] = E_k[i][k][1];
            y[k] = E_k[i][k][0];
        }
        printf("fiting");
        fit(k_points[i], x,y);
    }
    return 0;
}
