/********************************************************************************************************************************
Name: Lennard-Jones 2D
Author: Mahadi torabi and Md Shafayet Islam
Date: 2025-03-21
Description:
    A program to simulate an N number of particles in a 2D box with PBC using Lennard-Jones potential at zero pressure and 
    temperature. The program will use the LeapFrog algorithm to integrate the equations of motion
Usage:
    The program will be executed from the command line as follows:
    ./lj_3d

Output:
    trajectory.xyz: The positions of the particles. .xyz format can be visualized using VMD
    energy.dat: The program will also output the energy of the system.
    
Variables:
    The program will is driven by these following global variables:
    N - number of particles (default: 40)
    steps - number of steps (default: 1e6)
    dt - time step (default: 1e-4)
    L - box length (default: 10.0)
    rc - Lennard-Jonnes cut-off distance (default: 2.5)
    output_freq - output frequency (default: 200)
    lattice_type - Lattice type initialization('t' for triangular, 's' for square) (default: 't')


*********************************************************************************************************************************/ 


#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>


// Constants
double sigma = 1.0;
double epsilon = 1.0;
double mass = 1.0;
int N = 210;
int steps = 1e5;
double dt = 0.0001;
double L = 10.0;
double rc = 2.5;
int output_freq = 100;
char lattice_type = 'c'; // lattice type ('c' for cubic, 'h' for hexagonal)





// function to initialize the positions of the particles on a grid equally spaced within the box
void initialize_positions(int n, double L, double *x, double *y, double *z, char lattice_type)
{
    if (lattice_type == 'c')
    {
        int M = ceil(cbrt(n));
        double a = L / M;
        int k = 0;
        for (int i = 0; i < M; i++)
        {
            for (int j = 0; j < M; j++)
            {
                for (int l = 0; l < M; l++)
                {
                    if (k < n)
                    {
                        x[k] = (i + 0.5) * a;
                        y[k] = (j + 0.5) * a;
                        z[k] = (l + 0.5) * a;
                        k++;
                    }
                }
            }
        }
    }
    else if (lattice_type == 'h'){
       // Determine lattice constant 'a' by ensuring we fit enough atoms in the box
    double a = L / ceil(cbrt(n));  
    double c = sqrt(8.0 / 3.0) * a;  // HCP layer spacing

    // Compute the number of unit cells along each axis
    int nx = static_cast<int>(ceil(L / a));
    int ny = static_cast<int>(ceil(L / (a * sqrt(3) / 2.0)));
    int nz = static_cast<int>(ceil(L / c));

    int i = 0;
    for (int iz = 0; iz < nz && i < n; iz++) {
        for (int iy = 0; iy < ny && i < n; iy++) {
            for (int ix = 0; ix < nx && i < n; ix++) {
                if (i < n) {
                    // HCP stacking (A-B-A-B layers)
                    x[i] = ix * a + 0.5 * (iy % 2) * a;
                    y[i] = iy * a * sqrt(3) / 2.0;
                    z[i] = iz * c;

                    // Apply periodic boundary conditions
                    x[i] = x[i] - L * floor(x[i] / L);
                    y[i] = y[i] - L * floor(y[i] / L);
                    z[i] = z[i] - L * floor(z[i] / L);

                    // Check for NaN values
                    if (std::isnan(x[i]) || std::isnan(y[i]) || std::isnan(z[i])) {
                        std::cerr << "NaN detected at particle " << i << std::endl;
                        return;
                    }

                    i++;
                }
            }
        }
    }

    // Check if all particles were initialized
    if (i < n) {
        std::cerr << "Warning: Only " << i << " particles were initialized out of " << n << std::endl;
    }

    }
    
    else
    {
        std::cerr << "Invalid lattice type. Please enter 's' for simple cubic lattice or 'h' for hexagonal lattice." << std::endl;

    }



}





// function to initialize the velocities of the particles, with the center of mass velocity set to zero and the temperature set to 0
void initialize_velocities(int n, double *vx, double *vy, double *vz)
{
    double vxcm = 0.0;
    double vycm = 0.0;
    double vzcm = 0.0;
    for (int i = 0; i < n; i++)
    {
        //random velocity between -0.4 and 0.4
        vx[i] = -0.4 + 0.8 * rand() / RAND_MAX;
        vy[i] = -0.4 + 0.8 * rand() / RAND_MAX;
        vz[i] = -0.4 + 0.8 * rand() / RAND_MAX;


        vxcm += vx[i];
        vycm += vy[i];
        vzcm += vz[i];
    }
    vxcm /= n; // center of mass velocity
    vycm /= n;
    vzcm /= n;
    // set the center of mass velocity to zero
    for (int i = 0; i < n; i++)
    {
        vx[i] -= vxcm;
        vy[i] -= vycm;
        vz[i] -= vzcm;
    }
    
}


// function to calculate the distance between two particles with PBC
double distance(double x1, double y1, double z1, double x2, double y2, double z2, double L)
{
    double dx = x2 - x1;
    double dy = y2 - y1;
    double dz = z2 - z1;
    dx -= L * round(dx / L);
    dy -= L * round(dy / L);
    dz -= L * round(dz / L);
    return sqrt(dx * dx + dy * dy + dz * dz);
}

// function to calculate the distance between two particles without PBC
double distance_no_pbc(double x1, double y1, double z1, double x2, double y2, double z2)
{
    double dx = x2 - x1;
    double dy = y2 - y1;
    double dz = z2 - z1;
    return sqrt(dx * dx + dy * dy + dz * dz);
}


// function to calculate the Lennard-Jones potential
double lj_potential(double r)
{
    return 4.0 * epsilon * (pow(sigma / r, 12) - pow(sigma / r, 6));
}

// function to calculate the Lennard-Jones force
double lj_force(double r)
{
    return 24.0 * epsilon * (2.0 * pow(sigma / r, 13) - pow(sigma / r, 7)) / r;
}

// function to calculate the total energy of the system
double total_energy(int n, double L, double rc, double *x, double *y, double *z, double *vx, double *vy, double *vz)
{
    double potential = 0.0;
    double r = 0.0;
    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            r = distance(x[i], y[i], z[i], x[j], y[j], z[j], L);
            if (r < rc)
            {
                potential += lj_potential(r);
            }
        }
    }
    double kinetic = 0.0;
    for (int i = 0; i < n; i++)
    {
        kinetic += 0.5 * mass * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
    }
    return potential + kinetic;
}

// function to calculate the total energy of the system without PBC
double total_energy_no_pbc(int n, double *x, double *y, double *z)
{
    double potential = 0.0;
    double r = 0.0;
    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            r = distance_no_pbc(x[i], y[i], z[i], x[j], y[j], z[j]);
            potential += lj_potential(r);
        }
    }
    return potential;
}


// function to calculate the forces on the particles
void calculate_forces(int n, double L, double rc, double *x, double *y, double *z, double *fx, double *fy, double *fz)
{
    double dx, dy, dz, r, f = 0.0;

    for (int i = 0; i < n; i++)
    {
        fx[i] = 0.0;
        fy[i] = 0.0;
        fz[i] = 0.0;
    }
    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            dx = x[j] - x[i];
            dy = y[j] - y[i];
            dz = z[j] - z[i];
            dx -= L * round(dx / L);
            dy -= L * round(dy / L);
            dz -= L * round(dz / L);
            r = sqrt(dx * dx + dy * dy + dz * dz);
            if (r < rc)
            {
                f = lj_force(r);
                fx[i] += f * dx / r;
                fy[i] += f * dy / r;
                fz[i] += f * dz / r;
                fx[j] -= f * dx / r;
                fy[j] -= f * dy / r;
                fz[j] -= f * dz / r;
            }
        }
    }
}


// function to calculate the forces on the particles without PBC
void calculate_forces_no_pbc(int n, double *x, double *y, double *z, double *fx, double *fy, double *fz)
{
    double dx, dy, dz, r, f = 0.0;

    for (int i = 0; i < n; i++)
    {
        fx[i] = 0.0;
        fy[i] = 0.0;
        fz[i] = 0.0;
    }
    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            dx = x[j] - x[i];
            dy = y[j] - y[i];
            dz = z[j] - z[i];
            r = sqrt(dx * dx + dy * dy + dz * dz);
            if (r < rc)
            {
                f = lj_force(r);
                fx[i] += f * dx / r;
                fy[i] += f * dy / r;
                fz[i] += f * dz / r;
                fx[j] -= f * dx / r;
                fy[j] -= f * dy / r;
                fz[j] -= f * dz / r;
            }
        }
    }
}


//function to update the positions of the particles
void update_positions(int n, double L, double dt, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *fx, double *fy, double *fz)
{
    for (int i = 0; i < n; i++)
    {
        x[i] += vx[i] * dt + 0.5 * fx[i] * dt * dt / mass;
        y[i] += vy[i] * dt + 0.5 * fy[i] * dt * dt / mass;
        z[i] += vz[i] * dt + 0.5 * fz[i] * dt * dt / mass;
        x[i] -= L * round(x[i] / L);
        y[i] -= L * round(y[i] / L);
        z[i] -= L * round(z[i] / L);
    }
}


// function to update the positions of the particles without PBC
void update_positions_no_pbc(int n, double dt, double *x, double *y, double *z, double *fx, double *fy, double *fz)
{
    for (int i = 0; i < n; i++)
    {
        x[i] -=  fx[i] * dt ;
        y[i] -=  fy[i] * dt ;
        z[i] -=  fz[i] * dt ;
    }
}

// function to update the velocities of the particles
void update_velocities(int n, double dt, double *vx, double *vy, double *vz, double *fx, double *fy, double *fz, double *fx1, double *fy1, double *fz1)
{
    for (int i = 0; i < n; i++)
    {
        vx[i] += 0.5 * (fx[i] + fx1[i]) * dt / mass;
        vy[i] += 0.5 * (fy[i] + fy1[i]) * dt / mass;
        vz[i] += 0.5 * (fz[i] + fz1[i]) * dt / mass;
    }
}

// function to run the simulation
void run_simulation(int n, int steps, double dt, double L, double rc, double T, int output_freq, int seed, double *x, double *y, double *z, double *vx, double *vy, double *vz,char lattice_type ,std::string basename)
{
    srand(seed);
    double *fx = new double[n];
    double *fy = new double[n];
    double *fz = new double[n];
    double *fx1 = new double[n];
    double *fy1 = new double[n];
    double *fz1 = new double[n];

    std::cout << "Running simulation..." << std::endl;

    //set_parameters(n, steps, dt, L, rc, T, output_freq, seed);
    initialize_positions(n, L, x, y, z, lattice_type);
    initialize_velocities(n, vx, vy, vz);
    std::string energy_filename = basename + "_energy.dat";
    std::string trajectory_filename = basename + "_trajectory.xyz";
    std::ofstream energy_file(energy_filename);
    std::ofstream trajectory_file(trajectory_filename);
    trajectory_file << n << std::endl;

    calculate_forces(n, L, rc, x, y, z, fx, fy, fz);

    for (int step = 0; step < steps; step++)
    {
        
        for (int i = 0; i < n; i++)
        {
            fx1[i] = fx[i];
            fy1[i] = fy[i];
            fz1[i] = fz[i];
        }

        update_positions(n, L, dt, x, y, z, vx, vy, vz, fx, fy, fz);
        calculate_forces(n, L, rc, x, y, z, fx, fy, fz);
        update_velocities(n, dt, vx, vy, vz, fx, fy, fz, fx1, fy1, fz1);
        
        if (step % output_freq == 0)
        {
            double energy = total_energy(n, L, rc, x, y, z, vx, vy, vz);
            energy_file << step << " " << energy << std::endl;
            
            //trajectory_file << "Step " << step << std::endl;
            for (int i = 0; i < n; i++)
            {
                trajectory_file << "H " << x[i] << " " << y[i] << " " << z[i] << std::endl;
            }
        }
    }
    energy_file.close();
    trajectory_file.close();
    delete[] fx;
    delete[] fy;
    delete[] fz;
    delete[] fx1;
    delete[] fy1;
    delete[] fz1;
}

// function to perform energy minimization on the last configuration of the simulation
void minimize_energy(int n, double rc, double *x, double *y, double *z, double tol  ,std::string basename)
{
    double *fx = new double[n];
    double *fy = new double[n];
    double *fz = new double[n];
    double dt = 0.001;
    double max_force = 1.0;
    int step = 0;
    const int MAX_STEPS = 1e7; // Prevent infinite loops
    double f = 0.0;

    std::cout << "Minimizing energy..." << std::endl;
    std::string minimized_energy_filename = basename + "_min_energy.dat"; 
    std::ofstream minimized_energy_file(minimized_energy_filename);

    double prev_energy = total_energy_no_pbc(n, x, y, z);

    while (max_force > tol && step < MAX_STEPS)
    {
        calculate_forces_no_pbc(n, x, y, z, fx, fy, fz);
        
        // Update positions
        update_positions_no_pbc(n, dt, x, y, z, fx, fy, fz);
        
        // Recalculate forces
        calculate_forces_no_pbc(n, x, y, z, fx, fy, fz);

        // Compute max force
        max_force = 0.0;
        for (int i = 0; i < n; i++)
        {
            f = sqrt(fx[i] * fx[i] + fy[i] * fy[i] + fz[i] * fz[i]);
            if (f > max_force) {
                max_force = f;
            }
        }

        // Adaptive time step
        dt = std::min(0.001, tol / std::max(max_force, 1.0));

        double current_energy = total_energy_no_pbc(n, x, y, z);
        if (step % output_freq == 0)
        {
            minimized_energy_file << step << " " << current_energy << std::endl;
        }


        if (std::abs(current_energy - prev_energy) < 1e-7)
        {
           
            break; // Exit if energy is not changing
        }

        prev_energy = current_energy;

        step++;
        
    }

    // Save energy to file
   
    minimized_energy_file.close();

    std::string minimized_trajectory_filename = basename + "_min_config.xyz";
    // Save final configuration in XYZ format
    std::ofstream minimized_trajectory_file(minimized_trajectory_filename);
    minimized_trajectory_file << n << "\n\n";  // Add comment line for XYZ format
    for (int i = 0; i < n; i++)
    {
        minimized_trajectory_file << "H " << x[i] << " " << y[i] << " " << z[i] << std::endl;
    }
    minimized_trajectory_file.close();

    delete[] fx;
    delete[] fy;
    delete[] fz;

    std::cout << "Minimization completed in " << step << " steps." << std::endl;
}





int main()
{
    double *x = new double[N];
    double *y = new double[N];
    double *z = new double[N];
    double *vx = new double[N];
    double *vy = new double[N];
    double *vz = new double[N];
    
    std::string outdir = "run/"; //If you want to change the output directory, don't forget to change the output directory in the makefile in "test"
    std::string basename = outdir+"lj_3d";


    auto start = std::chrono::high_resolution_clock::now();
    run_simulation(N, steps, dt, L, rc, 1.0, output_freq, 12345, x, y, z, vx, vy, vz, lattice_type ,basename);
    // End timing
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
 
    std::cout << "Simulation completed in " << elapsed.count() << " seconds." << std::endl;
 
    double tol = 1e-4;
    minimize_energy(N, rc, x, y, z, tol, basename);

    delete[] x;
    delete[] y;
    delete[] z;
    delete[] vx;
    delete[] vy;
    delete[] vz;

    return 0;
}