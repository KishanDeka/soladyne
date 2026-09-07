import numpy as np
import matplotlib.pyplot as plt

file_name = 'viscosity.dat'
data = np.genfromtxt(file_name)

#azimuths = np.linspace(3.1415927, 0.0, 297)
#zeniths = np.linspace(3.828, 6.96, 297)

rfrac = np.array([(data[i][0])/6.96 for i in range(257)])
nuD = np.array([data[i][1] for i in range(257)])       

fig = plt.figure()
plt.plot(rfrac,nuD)
plt.xlabel(r"$r/R_s$")
plt.ylabel(r"viscosity $(10^{12}~cm^2/s)$")
plt.savefig('viscosity.jpg',dpi=500, bbox_inches='tight')

