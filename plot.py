from cProfile import label
from cmath import cos, pi
from telnetlib import GA
from turtle import color
import matplotlib.pyplot as plt
import numpy as np
from pyparsing import line
import math

band = open('bandstrc.dat', "r")
polynum = open('polynum.dat', "r")
gapfile = open('gap.dat',"r") 
E = []
K = []
poly_coeff = []
gaps= []


def PolyNoms(x, coeffs):
    y =0
    for i in range(len(coeffs)):
        y += coeffs[i] * x**i

    return y
def Cos(x,n,coeff):
    A = (coeff[0] + coeff[1])/2
    B = (coeff[0] - coeff[1])/2
    return A  * B * math.cos(math.pi * x)


for row in band:
    row = row.split('\t')
    try:
        E.append(float(row[0]))
        K.append(float(row[1]))

    except ValueError:
        continue
for row in gapfile:
    gaps.append([])
    row = row.split('\t')
    try:
        gaps[-1].append(float(row[0]))
        gaps[-1].append(float(row[1]))
    except ValueError:
        continue
print(gaps)
line_polynum = [line.split("\t") for line in polynum]
for line in line_polynum:
    poly_coeff.append([float(x) for x in line[:-1]])


k_points = np.linspace(-1, 1, 100)
Y_polynum= []
Y_cos=[]
for co_i in range(len(poly_coeff)):
    Y_polynum.append([])
    for k in k_points:
        Y_polynum[-1].append(PolyNoms(k,poly_coeff[co_i]))

for i in range(len(gaps)):
    Y_cos.append([])
    for k in k_points:
        Y_cos[-1].append( Cos(k,i,gaps[i]) )



for band in Y_polynum:
    plt.plot(k_points, band, 'b-',linewidth=0.8)
for band in Y_cos:
    plt.plot(k_points, band, 'g-',linewidth=0.8)


Gaps_polynum = []
Gaps_cos = []

for k in k_points:
    Gaps_polynum.append( PolyNoms(k,poly_coeff[1]) - PolyNoms(k,poly_coeff[0]) ) 
    Gaps_cos.append( Cos(k,1,gaps[1]) - Cos(k,2,gaps[0]) )
print(np.average(Gaps_polynum))
#print(sum(gap)/len(gap))
plt.scatter(K, E, color='r', s=1.5, label="a,b = 10, 6")
plt.xlabel(r'K-points[${\pi}$/a]', fontsize=12)
plt.ylabel('E [eV]', fontsize=12)
#print(Gaps_polynum)
text_kwargs1 = dict(ha='center', va='center', fontsize=12, color='C3')
text_kwargs2 = dict(ha='center', va='center', fontsize=12, color='C2')
textY0 = 1.4
textDx = -0.1
plt.text(-0.7, textY0, '{}< Gap <{}'.format(round(min(Gaps_polynum),2),round(max(Gaps_polynum),2)), **text_kwargs1)
textY0 += textDx
plt.text(-0.7, textY0, 'Average:{}'.format(round(np.average(Gaps_polynum),2)), **text_kwargs1)
textY0 += textDx
plt.text(-0.7, textY0, '{}< Gap <{}'.format(round(min(Gaps_cos),2),round(max(Gaps_cos),2)), **text_kwargs2)
textY0 += textDx
plt.text(-0.7, textY0, 'Average:{}'.format(round(np.average(Gaps_cos),2)), **text_kwargs2)

plt.title('KronigPenney', fontsize=20)
#plt.ylim(0.8,1.4)
plt.xticks(np.arange(-1.0, +1.01, .5))
#plt.yticks(np.arange(0.8, +1.7, .1))
plt.legend()
plt.savefig("band3" ,dpi=600)
# plt.show()
