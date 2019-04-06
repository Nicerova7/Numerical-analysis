#Taylor integration
import numpy as np
import matplotlib.pyplot as plt

def taylor(x):

    t = x*1.0
    z = x*1.0
    cont = 1
    
    while(True):
        t = -(x**2)*t*(2*cont-1)/(cont*(2*cont+1))
        cont += 1
        z += t
        if (np.abs(t) < 1E-8 ):
            break

    return z



val_x = np.linspace(0,1,10+1)
z = np.zeros(len(val_x))
for x in range(len(val_x)):
    z[x] = taylor(val_x[x])

plt.plot(val_x,z)
#print(z)
plt.show()
