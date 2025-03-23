/********************************************************************************************************************************
Name: Lennard-Jones 2D
Author: Mahadi torabi and Md Shafayet Islam
Date: 2025-03-21
Description:
    A program to simulate an N number of particles in a 2D box with PBC using Lennard-Jones potential at zero pressure and 
    temperature. The program will use the LeapFrog algorithm to integrate the equations of motion
Arguments:
    The program will take the following arguments:

Output:
    trajectory.xyz: The positions of the particles. .xyz format can be visualized using VMD
    energy.dat: The program will also output the energy of the system.
*/ 


#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

// Constants
constexpr double sigma = 1.0;
constexpr double epsilon = 1.0;
constexpr double mass = 1.0;
//constexpr double kB = 1.0;
constexpr int N = 20; // number of particles
constexpr int steps = 1e6; // number of steps
constexpr double dt = 1e-4; // time step
constexpr double L = 10.0; // box length
constexpr double rc = 2.5; // cut-off distance
//constexpr double T = 1.0;
constexpr int output_freq = 100; // output frequency


// function to initialize the positions of the particles in a lattice

constexpr void initialize_positions(int n, double L, double *x, double *y, char lattice_type)
{
    if (lattice_type == 't')
    {
        int nx = ceil(sqrt(n));
        int ny = (n + nx - 1) / nx; 
        double dx = L / nx;
        double dy = L / ny;
        int i = 0;
        for (int iy = 0; iy < ny; iy++)
        {
            for (int ix = 0; ix < nx; ix++)
            {
                if (i < n)
                {
                    x[i] = (ix + 0.5 * iy) * dx;
                    y[i] = iy * dy;
                    i++;
                }
            }
        }
    }
    else if (lattice_type == 's')
    {
        int nx = ceil(sqrt(n));
        int ny = (n + nx - 1) / nx; 
        double dx = L / nx;
        double dy = L / ny;
        int i = 0;
        for (int iy = 0; iy < ny; iy++)
        {
            for (int ix = 0; ix < nx; ix++)
            {
                if (i < n)
                {
                    x[i] = ix * dx;
                    y[i] = iy * dy;
                    i++;
                }
            }
        }
    }
    else
    {
        std::cerr << "Invalid lattice type. Please enter 't' for triangular lattice or 's' for square lattice." << std::endl;
    }
}



// function to initialize the velocities of the particles, with the center of mass velocity set to zero 
constexpr void initialize_velocities(int n, double *vx, double *vy)
{
    double vxcm = 0.0;
    double vycm = 0.0;
    for (int i = 0; i < n; i++)
    {
        //random velocity between -0.1 and 0.1
        vx[i] = 1e-1 * (2.0 * rand() / RAND_MAX - 1.0);
        vy[i] = 1e-1 * (2.0 * rand() / RAND_MAX - 1.0); // random velocity between -0.1 and 0.1

        vxcm += vx[i];
        vycm += vy[i];
    }
    vxcm /= n; // center of mass velocity
    vycm /= n;


    // set the center of mass velocity to zero
    for (int i = 0; i < n; i++)
    {
        vx[i] -= vxcm;
        vy[i] -= vycm;
    }

    
}


// function to calculate the distance between two particles with PBC
double distance(double x1, double y1, double x2, double y2, double L)
{
    double dx = x2 - x1;
    double dy = y2 - y1;
    dx -= L * round(dx / L);
    dy -= L * round(dy / L);
    return sqrt(dx * dx + dy * dy);
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
double total_energy(int n, double L, double rc, double *x, double *y, double *vx, double *vy)
{
    double energy = 0.0;
    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            double r = distance(x[i], y[i], x[j], y[j], L);
            if (r < rc)
            {
                energy += lj_potential(r);
            }
        }
    }
    double kinetic = 0.0;
    for (int i = 0; i < n; i++)
    {
        kinetic += 0.5 * mass * (vx[i] * vx[i] + vy[i] * vy[i]);
    }
    return energy + kinetic;
}



// function to calculate the forces on the particles
void calculate_forces(int n, double L, double rc, double *x, double *y, double *vx, double *vy, double *fx, double *fy, bool use_pbc)
{
    if (use_pbc==true)
    {
   
    
        for (int i = 0; i < n; i++)
        {
            fx[i] = 0.0;
            fy[i] = 0.0;
        }
        for (int i = 0; i < n; i++)
        {
            for (int j = i + 1; j < n; j++)
            {
                double dx = x[j] - x[i];
                double dy = y[j] - y[i];
                dx -= L * round(dx / L);
                dy -= L * round(dy / L);
                double r = sqrt(dx * dx + dy * dy);
                if (r < rc)
                {
                    double f = lj_force(r);
                    fx[i] += f * dx / r;
                    fy[i] += f * dy / r;
                    fx[j] -= f * dx / r;
                    fy[j] -= f * dy / r;
                }
            }
        }
    }
    else{
        for (int i = 0; i < n; i++)
        {
            fx[i] = 0.0;
            fy[i] = 0.0;
        }
        for (int i = 0; i < n; i++)
        {
            for (int j = i + 1; j < n; j++)
            {
                double dx = x[j] - x[i];
                double dy = y[j] - y[i];

                double r = sqrt(dx * dx + dy * dy);
                if (r < rc)
                {
                    double f = lj_force(r);
                    fx[i] += f * dx / r;
                    fy[i] += f * dy / r;
                    fx[j] -= f * dx / r;
                    fy[j] -= f * dy / r;
                }
            }
        }
    }

}



//function to update the positions of the particles

void update_positions(int n, double L, double dt, double *x, double *y, double *vx, double *vy, double *fx, double *fy)
{
    for (int i = 0; i < n; i++)
    {
        x[i] += vx[i] * dt + 0.5 * fx[i] * dt * dt / mass;
        y[i] += vy[i] * dt + 0.5 * fy[i] * dt * dt / mass;
        x[i] -= L * floor(x[i] / L);
        y[i] -= L * floor(y[i] / L);
    }
}

// function to update the velocities of the particles
void update_velocities(int n, double dt, double *vx, double *vy, double *fx, double *fy, double *fx1, double *fy1)
{
    for (int i = 0; i < n; i++)
    {
        vx[i] += 0.5 * (fx[i] + fx1[i]) * dt / mass;
        vy[i] += 0.5 * (fy[i] + fy1[i]) * dt / mass;
    }
}

// function to run the simulation
void run_simulation(int n, int steps, double dt, double L, double rc, double T, int output_freq, int seed, double *x, double *y, double *vx, double *vy,char lattice_type, std::string basename)
{
    srand(seed);
    double *fx = new double[n];
    double *fy = new double[n];
    double *fx1 = new double[n];
    double *fy1 = new double[n];
    //set_parameters(n, steps, dt, L, rc, T, output_freq, seed);
    initialize_positions(n, L, x, y, lattice_type);
    initialize_velocities(n, vx, vy);
    std::string energy_filename = basename + "_energy.dat";
    std::string trajectory_filename = basename + "_trajectory.xyz";
    std::ofstream energy_file(energy_filename);
    std::ofstream trajectory_file(trajectory_filename);
    trajectory_file << n << std::endl;
    for (int step = 0; step < steps; step++)
    {
        calculate_forces(n, L, rc, x, y, vx, vy, fx, fy, true);
        for (int i = 0; i < n; i++)
        {
            fx1[i] = fx[i];
            fy1[i] = fy[i];
        }
        update_positions(n, L, dt, x, y, vx, vy, fx, fy);
        calculate_forces(n, L, rc, x, y, vx, vy, fx, fy, true);
        update_velocities(n, dt, vx, vy, fx, fy, fx1, fy1);
        
        if (step % output_freq == 0)
        {
            double energy = total_energy(n, L, rc, x, y, vx, vy);
            energy_file << step << " " << energy << std::endl;
            
            //trajectory_file << "Step " << step << std::endl;
            for (int i = 0; i < n; i++)
            {
                trajectory_file << "H " << x[i] << " " << y[i] << " 0.0" << std::endl;
            }
        }
    }
    energy_file.close();
    trajectory_file.close();
    delete[] fx;
    delete[] fy;
    delete[] fx1;
    delete[] fy1;
}


void minimize_energy(int n, double L, double rc, double *x, double *y, double *vx, double *vy, bool use_pbc, double tol, std::string base)
{
    double *fx = new double[n];
    double *fy = new double[n];
    double dt = 0.01;
    double max_force = 1.0;
    int step = 0;
    const int MAX_STEPS = 1e4; 
    std::string minimized_energy_filename = base + "_min_energy.dat";
    std::ofstream minimized_energy_file(minimized_energy_filename);

    while (max_force > tol && step < MAX_STEPS)
    {
        calculate_forces(n, L, rc, x, y, vx, vy, fx, fy, use_pbc);
        
        // Update positions
        //update_positions(n, L, dt, x, y, vx, vy, fx, fy);
        // Manually update positions
        for (int i = 0; i < n; i++)
        {
            x[i] -= dt * fx[i];
            y[i] -= dt * fy[i];

        }
        
        // Recalculate forces
        calculate_forces(n, L, rc, x, y, vx, vy, fx, fy, use_pbc);

        // Compute min force
        max_force = 0.0;
        for (int i = 0; i < n; i++)
        {
            double f = sqrt(fx[i] * fx[i] + fy[i] * fy[i]);
            if (f > max_force) {
                max_force = f;
            }
        }

        // Adaptive time step
        dt = std::min(0.001, tol / std::max(max_force, 1.0));
        // Save energy to file
        minimized_energy_file << step << " " << total_energy(n, L, rc, x, y, vx, vy) << std::endl;

        step++;
    }

    minimized_energy_file.close();

    // Save final configuration in XYZ format
    std::string minimized_trajectory_filename = base + "_min_trajectory.xyz";
    std::ofstream minimized_trajectory_file(minimized_trajectory_filename);
    minimized_trajectory_file << n << "\n\n";  // Add comment line for XYZ format
    for (int i = 0; i < n; i++)
    {
        minimized_trajectory_file << "H " << x[i] << " " << y[i] << " 0.0" << std::endl;
    }
    minimized_trajectory_file.close();

    delete[] fx;
    delete[] fy;

}





int main()
{
    double *x = new double[N];
    double *y = new double[N];
    double *vx = new double[N];
    double *vy = new double[N];
    char lattice_type = 't';
    std::string outdir = "run/"; //If you want to change the output directory, don't forget to change the output directory in the makefile in "test"
    std::string basename = outdir+"lj_2d";
    run_simulation(N, steps, dt, L, rc, 1.0, output_freq, 12345, x, y, vx, vy, lattice_type, basename);
    
    bool use_pbc = false;
    minimize_energy(N, L, rc, x, y, vx, vy, use_pbc, 1e-1, basename );

    delete[] x;
    delete[] y;
    delete[] vx;
    delete[] vy;

    return 0;
}