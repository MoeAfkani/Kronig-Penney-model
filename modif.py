from cmath import sqrt
import matplotlib.pyplot as plt
import numpy as np
import math as m

def cosK(alpha):
    y = []
    beta = sqrt(1-alpha**2)
    for i in np.linspace(0,alpha,10):
        for j in np.linspace(0,beta,10):
           y.append( m.cos(j)*m.cosh(i) - m.sin(j) * m.sinh(i) / (2*alpha*beta) )


    return y

X = np.linspace(0,1,100)

for i in range (100):
    plt.plot(X , cosK(2.5))
plt.show()

