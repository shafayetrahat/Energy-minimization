/********************************************************************************************************************************
Name: Lennard-Jones 2D
Author: Mahadi torabi and Md Shafayet Islam
Date: 2025-03-21
Description:
    A program to simulate an N number of particles in a 2D box with PBC using Lennard-Jones potential at zero pressure and 
    temperature. The program will use the LeapFrog algorithm to integrate the equations of motion
Usage:
    The program will be executed from the command line as follows:
    ./lj_2d_variable_box

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
#include <string>
#include <chrono> // Include the chrono library for timing


// Constants
double sigma = 1.0;
double epsilon = 1.0;
double mass = 1.0;
int steps = 1e6; // number of steps
double dt = 1e-4; // time step
double rc = 2.5; // cut-off distance
int output_freq = 200; // output frequency
char lattice_type = 's'; // lattice type ('t' for triangular, 's' for square)


// function to initialize the positions of the particles in a lattice
void initialize_positions(int n, double L, double *x, double *y, char lattice_type)
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

void initialize_velocities(int n, double *vx, double *vy)
{
    double vxcm = 0.0;
    double vycm = 0.0;
    for (int i = 0; i < n; i++)
    {
        //random velocity between -0.1 and 0.1
        vx[i] = -0.1 + 0.2 * rand() / RAND_MAX;
        vy[i] = -0.1 + 0.2 * rand() / RAND_MAX;
        

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
    double potential = 0.0;
    double r = 0.0;

    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            r = distance(x[i], y[i], x[j], y[j], L);
            if (r < rc)
            {
                potential += lj_potential(r);
            }

        }

    }

    double kinetic = 0.0;
    for (int i = 0; i < n; i++)
    {
        kinetic += 0.5 * mass * (vx[i] * vx[i] + vy[i] * vy[i]);
    }
    return potential + kinetic;
}



// function to calculate the forces on the particles
void calculate_forces(int n, double L, double rc, double *x, double *y, double *vx, double *vy, double *fx, double *fy)
{

    double dx, dy, r, f = 0.0;

    for (int i = 0; i < n; i++)
    {
        fx[i] = 0.0;
        fy[i] = 0.0;
    }
    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            dx = x[j] - x[i];
            dy = y[j] - y[i];
            dx -= L * round(dx / L);
            dy -= L * round(dy / L);
            r = sqrt(dx * dx + dy * dy);
            if (r < rc)
            {
                f = lj_force(r);
                fx[i] += f * dx / r;
                fy[i] += f * dy / r;
                fx[j] -= f * dx / r;
                fy[j] -= f * dy / r;
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
        x[i] -= L * round(x[i] / L);
        y[i] -= L * round(y[i] / L);
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
void run_simulation(int n, int steps, double dt, double L, double rc, double T, int output_freq, int seed, double *x, double *y, double *vx, double *vy, char lattice_type)
{
    srand(seed);
    double *fx = new double[n];
    double *fy = new double[n];
    double *fx1 = new double[n];
    double *fy1 = new double[n];

    
    initialize_positions(n, L, x, y, lattice_type);
    initialize_velocities(n, vx, vy);
    

    calculate_forces(n, L, rc, x, y, vx, vy, fx, fy);

    for (int step = 0; step < steps; step++)
    {
        
        for (int i = 0; i < n; i++)
        {
            fx1[i] = fx[i];
            fy1[i] = fy[i];
        }

        update_positions(n, L, dt, x, y, vx, vy, fx, fy);
        calculate_forces(n, L, rc, x, y, vx, vy, fx, fy);
        update_velocities(n, dt, vx, vy, fx, fy, fx1, fy1);
        
       
    }
    
    delete[] fx;
    delete[] fy;
    delete[] fx1;
    delete[] fy1;
}








int main()
{
    double L = 10.0; // box length
    int N = 30; // number of particles
    double density = N / (L * L); // density


    double *x = new double[N];
    double *y = new double[N];
    double *vx = new double[N];
    double *vy = new double[N];

    




    //open a file to save the time for each box length
    std::ofstream file("run/time_vs_box_2d.dat");

    for(int i = 0; i < 18; i++){

       auto start = std::chrono::high_resolution_clock::now();

        run_simulation(N, steps, dt, L, rc, 1.0, output_freq, 12345, x, y, vx, vy, lattice_type);
  
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        file << L << " " << N << " "<< elapsed.count() << std::endl;
        L += 2.0;
        N = static_cast<int>(density * L * L);
        delete[] x;
        delete[] y;
        delete[] vx;
        delete[] vy;
        x = new double[N];
        y = new double[N];
        vx = new double[N];
        vy = new double[N];
    
    }

    file.close();

    
   
    delete[] x;
    delete[] y;
    delete[] vx;
    delete[] vy;

    return 0;
}