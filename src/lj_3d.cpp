#include <iostream>
#include <fstream>
#include <cmath>





// A program to simulate an N number of particles in a 2D box with PBC using Lennard-Jones potential at zero pressure and temperature
// The program will use the LeapFrog algorithm to integrate the equations of motion
// The program will output the positions of the particles in a file called "trajectory.xyz" which can be visualized using VMD
// The program will also output the energy of the system in a file called "energy.dat"






// Constants
constexpr double sigma = 1.0;
constexpr double epsilon = 1.0;
constexpr double mass = 1.0;
//constexpr double kB = 1.0;
constexpr int N = 60;
constexpr int steps = 1e5;
constexpr double dt = 0.001;
constexpr double L = 10.0;
constexpr double rc = 2.5;
//constexpr double T = 1.0;
constexpr int output_freq = 20;


// function to set the parameters of the simulation

// void set_parameters(int &n, int &steps, double &dt, double &L, double &rc, double &T, int &output_freq, int &seed)
// {
//     std::cout << "Enter the number of particles: ";
//     std::cin >> n;
//     std::cout << "Enter the number of steps: ";
//     std::cin >> steps;
//     std::cout << "Enter the time step: ";
//     std::cin >> dt;
//     std::cout << "Enter the box length: ";
//     std::cin >> L;
//     std::cout << "Enter the cut-off distance: ";
//     std::cin >> rc;
//     std::cout << "Enter the temperature: ";
//     std::cin >> T;
//     std::cout << "Enter the output frequency: ";
//     std::cin >> output_freq;
//     std::cout << "Enter the seed for random number generator: ";
//     std::cin >> seed;   
        
// }


// function to initialize the positions of the particles on a grid equally spaced within the box
constexpr void initialize_positions(int n, double L, double *x, double *y, double *z)
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


// function to initialize the positions of the particles in a triangular lattice within the box
// constexpr void initialize_positions(int n, double L, double *x, double *y, double *z)
// {
//     int M = ceil(sqrt(n));
//     double a = L / M;
//     int k = 0;
//     for (int i = 0; i < M; i++)
//     {
//         for (int j = 0; j < M; j++)
//         {
//             if (k < n)
//             {
//                 x[k] = (i + 0.5) * a;
//                 y[k] = (j + 0.5) * a;
//                 z[k] = 0.5 * a;
//                 k++;
//             }
//         }
//     }
// } 



// function to initialize the velocities of the particles, with the center of mass velocity set to zero and the temperature set to 0
constexpr void initialize_velocities(int n, double *vx, double *vy, double *vz)
{
    double vxcm = 0.0;
    double vycm = 0.0;
    double vzcm = 0.0;
    for (int i = 0; i < n; i++)
    {
        //random velocity between 0.5 and -0.5
        vx[i] = 0.5 - rand() / (RAND_MAX + 1.0);
        vy[i] = 0.5 - rand() / (RAND_MAX + 1.0);
        vz[i] = 0.5 - rand() / (RAND_MAX + 1.0);

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
    double energy = 0.0;
    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            double r = distance(x[i], y[i], z[i], x[j], y[j], z[j], L);
            if (r < rc)
            {
                energy += lj_potential(r);
            }
        }
    }
    double kinetic = 0.0;
    for (int i = 0; i < n; i++)
    {
        kinetic += 0.5 * mass * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
    }
    return energy + kinetic;
}


// function to calculate the forces on the particles
void calculate_forces(int n, double L, double rc, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *fx, double *fy, double *fz, bool use_pbc)
{
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
            double dx = x[j] - x[i];
            double dy = y[j] - y[i];
            double dz = z[j] - z[i];
            if (use_pbc)
            {
            dx -= L * round(dx / L);
            dy -= L * round(dy / L);
            dz -= L * round(dz / L);
            }
            double r = sqrt(dx * dx + dy * dy + dz * dz);
            if (r < rc)
            {
                double f = lj_force(r);
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
        x[i] -= L * floor(x[i] / L);
        y[i] -= L * floor(y[i] / L);
        z[i] -= L * floor(z[i] / L);
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
void run_simulation(int n, int steps, double dt, double L, double rc, double T, int output_freq, int seed, double *x, double *y, double *z, double *vx, double *vy, double *vz, bool use_pbc, std::string basename)
{
    srand(seed);
    double *fx = new double[n];
    double *fy = new double[n];
    double *fz = new double[n];
    double *fx1 = new double[n];
    double *fy1 = new double[n];
    double *fz1 = new double[n];
    //set_parameters(n, steps, dt, L, rc, T, output_freq, seed);
    initialize_positions(n, L, x, y, z);
    initialize_velocities(n, vx, vy, vz);
    std::string energy_filename = basename + "_energy.dat";
    std::string trajectory_filename = basename + "_trajectory.xyz";
    std::ofstream energy_file(energy_filename);
    std::ofstream trajectory_file(trajectory_filename);
    trajectory_file << n << std::endl;
    for (int step = 0; step < steps; step++)
    {
        calculate_forces(n, L, rc, x, y, z, vx, vy, vz, fx, fy, fz, use_pbc);
        for (int i = 0; i < n; i++)
        {
            fx1[i] = fx[i];
            fy1[i] = fy[i];
            fz1[i] = fz[i];
        }
        update_positions(n, L, dt, x, y, z, vx, vy, vz, fx, fy, fz);
        calculate_forces(n, L, rc, x, y, z, vx, vy, vz, fx, fy, fz, use_pbc);
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
void minimize_energy(int n, double L, double rc, double *x, double *y, double *z, double *vx, double *vy, double *vz, bool use_pbc, std::string basename)
{
    double *fx = new double[n];
    double *fy = new double[n];
    double *fz = new double[n];
    double dt = 0.01;
    double tol = 1.0e-3;
    double max_force = 1.0;
    int step = 0;
    const int MAX_STEPS = 1e7; // Prevent infinite loops
    std::string minimized_energy_filename = basename + "_min_energy.dat";
    
    std::ofstream minimized_energy_file(minimized_energy_filename);
    while (max_force > tol && step < MAX_STEPS)
    {
        calculate_forces(n, L, rc, x, y, z, vx, vy, vz, fx, fy, fz, use_pbc);
        
        // Manually update positions
        for (int i = 0; i < n; i++)
        {
            x[i] -= dt * fx[i];
            y[i] -= dt * fy[i];
            z[i] -= dt * fz[i];
        }
        
        // Recalculate forces
        calculate_forces(n, L, rc, x, y, z, vx, vy, vz, fx, fy, fz, use_pbc);

        // Compute max force
        max_force = 0.0;
        for (int i = 0; i < n; i++)
        {
            double f = sqrt(fx[i] * fx[i] + fy[i] * fy[i] + fz[i] * fz[i]);
            if (f > max_force) {
                max_force = f;
            }
        }

        // Adaptive time step
        dt = std::min(0.001, tol / std::max(max_force, 1.0));
        minimized_energy_file << step << " " << total_energy(n, L, rc, x, y, z, vx, vy, vz) << std::endl;
        step++;
    }

    // Save energy to file
   
    minimized_energy_file.close();
    std::string minimized_trajectory_filename = basename + "_min_trajectory.xyz";
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
}





int main()
{
    double *x = new double[N];
    double *y = new double[N];
    double *z = new double[N];
    double *vx = new double[N];
    double *vy = new double[N];
    double *vz = new double[N];
    bool use_pbc = true;
    std::string outdir = "run/"; //If you want to change the output directory, don't forget to change the output directory in the makefile in "test"
    std::string basename = outdir+"lj_3d";

    run_simulation(N, steps, dt, L, rc, 1.0, output_freq, 12345, x, y, z, vx, vy, vz, use_pbc,basename);
    use_pbc = false;
    minimize_energy(N, L, rc, x, y, z, vx, vy, vz,use_pbc,basename);

    delete[] x;
    delete[] y;
    delete[] z;
    delete[] vx;
    delete[] vy;
    delete[] vz;

    return 0;
}