from telnetlib import GA
from turtle import color
import matplotlib.pyplot as plt
import numpy as np
from pyparsing import line

band = open('bandstrc.dat', "r")
polynum = open('polynum.dat', "r")
E = []
K = []
poly_coeff = []


def PolyNoms(x, coeffs):
    y =0
    for i in range(len(coeffs)):
        y += coeffs[i] * x**i

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


k_points = np.linspace(-1, 1, 100)
bands= []
for co_i in range(len(poly_coeff)):
    bands.append([])
    for k in k_points:
        bands[-1].append(PolyNoms(k,poly_coeff[co_i]))

for band in bands:
    plt.plot(k_points, band, 'b-',linewidth=0.8)


def gap(k, m=0,n=1):
    return PolyNoms(k,poly_coeff[n]) -PolyNoms(k,poly_coeff[m])

Gaps = []
for k in k_points:
    Gaps.append( gap(k) )
print(np.average(Gaps))
#print(sum(gap)/len(gap))
plt.scatter(K, E, color='r', s=1.5, label="a,b = 8.9, 8")
plt.xlabel(r'K-points[${\pi}$/a]', fontsize=12)
plt.ylabel('E [eV]', fontsize=12)
#print(Gaps)
text_kwargs = dict(ha='center', va='center', fontsize=12, color='C1')
plt.text(0.0, 0.8, '{}< Gap <{}'.format(round(min(Gaps),2),round(max(Gaps),2)), **text_kwargs)
plt.text(0.0, 0.65, 'Average:{}'.format(round(np.average(Gaps),2)), **text_kwargs)

plt.title('KronigPenney', fontsize=20)
plt.xticks(np.arange(-1.0, +1.01, .5))
plt.yticks(np.arange(0.3, +1.5, .1))
plt.legend()
plt.savefig("band.png")
# plt.show()
