import numpy as np
import matplotlib.pyplot as plt
from scipy.special import*

x = np.arange(-20,5,0.1)

Z = airy(x)

Ai = Z[0]
Aip = Z[1]
Bi = Z[2]
Bip = Z[3]
f = Ai+Bi
fp = Aip+Bip

plt.plot(x,f)
plt.plot(x,fp)
axes = plt.gca()
axes.set_ylim([-5,5])
plt.show()

#You're all set for root finding
