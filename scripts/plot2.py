import matplotlib.pyplot as plt
import numpy as np

#num_particles = 30
data = np.loadtxt('time_vs_box_3d_linked_cell.dat')

x1 = data[:,0]
y1 = data[:,2]
n1 = data[:,1]



x1 = x1*x1*x1


# figure size
plt.figure(figsize=(10, 6))
plt.plot(x1, y1, "r", label= "Simulation Time") # total energy
plt.plot(x1,n1/32,"b", label="Number of Particles * (1/32)")
plt.xlabel('Box Size (Volume) ')

plt.ylabel('Simulation Time (sec)')
plt.legend()
plt.title('Simulation Time vs Box Size (Density=0.15)')
plt.savefig('time_vs_box_3d_linked_cell.png',dpi=300)


plt.show()
