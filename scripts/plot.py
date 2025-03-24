import matplotlib.pyplot as plt
import numpy as np

num_particles = 20
data = np.loadtxt('lj_3d_min_energy.dat')
x = data[:,0]
y = data[:,1]

# figure size
plt.figure(figsize=(10, 6))
plt.plot(x, y, "r", label= "Total energy") # total energy
#plt.plot(x, y/20, label= "Average energy per particle") # average energy per particle
plt.xlabel('Step')
#plt.xlim(0, 4e5)
# make sure the last one appears
#plt.xticks(np.arange(0, 1e6+1, 1e5))
#plt.xticks(rotation=45)
plt.ylabel('Energy')
plt.legend()
plt.title('Energy vs Step (Minimization)')
plt.savefig('min_energy_3d.png')
plt.show()

