import matplotlib.pyplot as plt
import numpy as np
from pyparsing import line

band = open('bandstrc.dat', "r")
polynum = open('polynum.dat', "r")
E = []
K = []
poly_coeff = []


def PolyCoefficients(x, coeffs):
    print('# This is a polynomial of order:', len(coeffs))
    y = []
    for x0 in range(len(x)):
        y.append(0)
        for i in range(len(coeffs)):
            y[x0] += coeffs[i] * x[x0]**i

    return y


for row in band:
    row = row.split('\t')
    try:
        E.append(float(row[0]))
        K.append(float(row[1]))

    except ValueError:
        continue

line_polynum = [line.split("\t") for line in polynum]
for line in line_polynum:
    poly_coeff.append([float(x) for x in line[:-1]])
print(poly_coeff)


x = np.linspace(-1, 1, 100)
for co_i in poly_coeff:
    plt.plot(x, PolyCoefficients(x, co_i), label="PolyFit")

plt.scatter(K, E, color='r', s=2, label='a,b = 8.9, 8')
plt.xlabel(r'K-points[${\pi}$/a]', fontsize=12)
plt.ylabel('E [eV]', fontsize=12)

plt.title('KronigPenney', fontsize=20)
plt.xticks(np.arange(-1.0, +1.01, .5))
plt.legend()
plt.savefig("band.png")
# plt.show()
