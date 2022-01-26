import matplotlib.pyplot as plt

band  = open('bandstrc.dat', "r")
E=[]
K=[]
for row in band:
    row = row.split('\t')
    E.append(float(row[0]))
    K.append(float(row[1]))
  
plt.scatter(K, E, color = 'r', s=0.8, label = 'a,b = 10, 8')
  
plt.xlabel('K-points', fontsize = 12)
plt.ylabel('E [eV]', fontsize = 12)
  
plt.title('KronigPenney', fontsize = 20)
plt.legend()
plt.savefig("band.png")
#plt.show()


print(band)