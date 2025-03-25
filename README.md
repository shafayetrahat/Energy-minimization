# Energy-minimization
## Description
## How to run
Command to compile and run the simulation: 
```
make clean
make
make test
```
Descriptions of the commands:
```
make clean
  - This will clean all of your binary files.
make cleanall
  - This will clean all of your binary files and output the results in the run directory.
make
  - This will create all the binaries in the bin/  directories.
make test
  - Will run all the simulations together.
```
### Simulation Parameters
```
N - number of particles ( default : 40)
steps - number of steps ( default : 1 e6 )
dt - time step ( default : 1e -4)
L - box length ( default : 10.0)
rc - Lennard - Jonnes cut - off distance ( default : 2.5)
output_freq - output frequency ( default : 200)
lattice_type - Lattice type initialization ( ’t ’ for triangular , ’s ’ for square ) ( default : ’t ’)
cell_size - per cell size
num_cells - number of cell
total_cells - total number of cells
```
### Output
The output will be saved in the run/ folder. However, if you enter the bin and want to run the individual binary, it will save the output in the child run/ directory (inside bin/ directory)
```
Naive approach output:
    lj_<dim>_trajectory.xyz: this file contains the trajectory, which can be visualized.
    lj_<dim>_energy.dat: this is for outputting the energy for each trajectory step for plotting.
    lj_<dim>_min_config.xyz: This file contains the configuration after the minimization has finished.
    lj_<dim>_min_energy.dat: This file contains the energy of the minimized system.
Linked-Cell algorithm outputs:
    lj_cell_<dim>_trajectory.xyz: this file contains the trajectory, which can be visualized.
    lj_cell_<dim>_energy.dat: this is for outputting the energy for each trajectory step for plotting.
    lj_cell_<dim>_min_config.xyz: This file contains the configuration after the minimization has finished.
    lj_cell_<dim>_min_energy.dat: This file contains the energy of the minimized system.
```
If you want to run a single binary code then you have to go to the "bin/" directory, create a "run/" directory there, and then run any binary file you want.
## Simulation Gallery
### Simulation of Lennard-Jones naive algorithm in 2d. The total particles are 20, 30, 40.

https://github.com/user-attachments/assets/f04f9054-8472-4362-8471-31572d6f0beb 

https://github.com/user-attachments/assets/9d4f7efb-3027-4f0a-a12a-9da4056e424d

https://github.com/user-attachments/assets/610a6f66-6751-4d8e-a949-b986ff2cff3f




### Simulation of Lennard-Jones linked cell-based algorithm in 2d. The total particles are 40.

https://github.com/user-attachments/assets/f7043450-dadd-4b29-921a-d30f72c11e4f


### Simulation of Lennard-Jones linked cell-based algorithm in 3d. The total particles are 210.

https://github.com/user-attachments/assets/872fc8ed-174d-4abe-b4a3-08ec52426adf

