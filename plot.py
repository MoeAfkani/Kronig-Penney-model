import matplotlib.pyplot as plt
import numpy as np

band  = open('bandstrc.dat', "r")
E=[]
K=[]
for row in band:
    row = row.split('\t')
    try:
        E.append(float(row[0]))
        K.append(float(row[1]))
        E.append(float(row[0]))
        K.append(-float(row[1]))
    except ValueError: continue
  
plt.scatter(K, E, color = 'r', s=0.8, label = 'a,b = 8.9, 8')
  
plt.xlabel(r'K-points[${\pi}$/a]', fontsize = 12)
plt.ylabel('E [eV]', fontsize = 12)
  
plt.title('KronigPenney', fontsize = 20)
plt.xticks(np.arange(-1, +1.01, .5))
plt.legend()
plt.savefig("band.png")
#plt.show()


print(band)