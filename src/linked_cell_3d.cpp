#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <vector>
#include <tuple>


// Constants
double sigma = 1.0;
double epsilon = 1.0;
double mass = 1.0;
int N = 210; // number of particles
int steps = 1e6; // number of steps
double dt = 1e-4; // time step
double L = 10.0; // box length
double rc = 2.4; // cut-off distance
int output_freq = 200; // output frequency

// Linked Cell Algorithm parameters
double cell_size = rc + 0.1;
int num_cells = static_cast<int>(L / cell_size);
int total_cells = num_cells * num_cells * num_cells;

struct Particle {
    double x, y, z, vx, vy, vz;
};

// Neighboring cell offsets
const std::vector<std::tuple<int, int, int>> neighbor_offsets = {
    std::make_tuple(-1, -1, -1), std::make_tuple(-1, -1, 0), std::make_tuple(-1, -1, 1),
    std::make_tuple(-1, 0, -1), std::make_tuple(-1, 0, 0), std::make_tuple(-1, 0, 1),
    std::make_tuple(-1, 1, -1), std::make_tuple(-1, 1, 0), std::make_tuple(-1, 1, 1),
    std::make_tuple(0, -1, -1), std::make_tuple(0, -1, 0), std::make_tuple(0, -1, 1),
    std::make_tuple(0, 0, -1), std::make_tuple(0, 0, 0), std::make_tuple(0, 0, 1),
    std::make_tuple(0, 1, -1), std::make_tuple(0, 1, 0), std::make_tuple(0, 1, 1),
    std::make_tuple(1, -1, -1), std::make_tuple(1, -1, 0), std::make_tuple(1, -1, 1),
    std::make_tuple(1, 0, -1), std::make_tuple(1, 0, 0), std::make_tuple(1, 0, 1),
    std::make_tuple(1, 1, -1), std::make_tuple(1, 1, 0), std::make_tuple(1, 1, 1)
};

// Function to initialize positions
void initialize_positions(int n, double L, Particle* particles) {
    int nx = static_cast<int>(ceil(pow(static_cast<double>(n), 1.0 / 3.0)));
    int ny = nx;
    int nz = nx;
    double dx = L / nx;
    double dy = L / ny;
    double dz = L / nz;
    int i = 0;
    for (int iz = 0; iz < nz; iz++) {
        for (int iy = 0; iy < ny; iy++) {
            for (int ix = 0; ix < nx; ix++) {
                if (i < n) {
                    // Add small random displacement (0.5% of cell size)
                    double rand_x = (2.0 * rand() / RAND_MAX - 1.0) * 0.005 * dx;
                    double rand_y = (2.0 * rand() / RAND_MAX - 1.0) * 0.005 * dy;
                    double rand_z = (2.0 * rand() / RAND_MAX - 1.0) * 0.005 * dz;
                    
                    particles[i].x = (ix + 0.5) * dx + rand_x;
                    particles[i].y = (iy + 0.5) * dy + rand_y;
                    particles[i].z = (iz + 0.5) * dz + rand_z;
                    i++;
                }
            }
        }
    }
}

// Function to initialize velocities
void initialize_velocities(int n, Particle* particles) {
    double vxcm = 0.0, vycm = 0.0, vzcm = 0.0;
    for (int i = 0; i < n; i++) {

        // Random velocities between -0.4 and 0.4
        // particles[i].vx = -0.4 + 0.8 * rand() / RAND_MAX;
        // particles[i].vy = -0.4 + 0.8 * rand() / RAND_MAX;
        // particles[i].vz = -0.4 + 0.8 * rand() / RAND_MAX;

        // Random velocities between -0.1 and 0.1
        particles[i].vx = -0.1 + 0.2 * rand() / RAND_MAX;
        particles[i].vy = -0.1 + 0.2 * rand() / RAND_MAX;
        particles[i].vz = -0.1 + 0.2 * rand() / RAND_MAX;


        vxcm += particles[i].vx;
        vycm += particles[i].vy;
        vzcm += particles[i].vz;
    }
    vxcm /= n;
    vycm /= n;
    vzcm /= n;
    for (int i = 0; i < n; i++) {
        particles[i].vx -= vxcm;
        particles[i].vy -= vycm;
        particles[i].vz -= vzcm;
    }
}

// Function to calculate distance with PBC
double distance(double x1, double y1, double z1, double x2, double y2, double z2, double L) {
    double dx = x2 - x1;
    double dy = y2 - y1;
    double dz = z2 - z1;
    dx -= L * round(dx / L);
    dy -= L * round(dy / L);
    dz -= L * round(dz / L);
    return sqrt(dx * dx + dy * dy + dz * dz);
}

// Lennard-Jones potential
double lj_potential(double r) {
    return 4.0 * epsilon * (pow(sigma / r, 12) - pow(sigma / r, 6));
}

// Lennard-Jones force
double lj_force(double r) {
    return 24.0 * epsilon * (2.0 * pow(sigma / r, 13) - pow(sigma / r, 7)) / r;
}

// Function to build the cell list
void build_cell_list(int n, double L, Particle* particles, std::vector<std::vector<int>>& cell_list) {
    for (auto& cell : cell_list) cell.clear();

    for (int i = 0; i < n; i++) {
        // Apply periodic boundary conditions using fmod
        particles[i].x = fmod(particles[i].x + L, L);
        particles[i].y = fmod(particles[i].y + L, L);
        particles[i].z = fmod(particles[i].z + L, L);

        // Clamp x, y, z to avoid indexing out of range
        particles[i].x = std::min(particles[i].x, L - 1e-12);
        particles[i].y = std::min(particles[i].y, L - 1e-12);
        particles[i].z = std::min(particles[i].z, L - 1e-12);

        // Calculate cell indices
        int cell_x = static_cast<int>(particles[i].x / cell_size);
        int cell_y = static_cast<int>(particles[i].y / cell_size);
        int cell_z = static_cast<int>(particles[i].z / cell_size);
        if(cell_x >= num_cells) cell_x = num_cells - 1;
        if(cell_y >= num_cells) cell_y = num_cells - 1;
        if(cell_z >= num_cells) cell_z = num_cells - 1;

        // Calculate cell index
        int cell_index = cell_x + cell_y * num_cells + cell_z * num_cells * num_cells;
        cell_list[cell_index].push_back(i);
    }
}

// Function to calculate forces using linked cell algorithm
void calculate_forces_linked_cell(int n, double L ,Particle* particles, double* fx, double* fy, double* fz, const std::vector<std::vector<int>>& cell_list) {
    std::fill(fx, fx + n, 0.0);
    std::fill(fy, fy + n, 0.0);
    std::fill(fz, fz + n, 0.0);

    for (int cell_x = 0; cell_x < num_cells; cell_x++) {
        for (int cell_y = 0; cell_y < num_cells; cell_y++) {
            for (int cell_z = 0; cell_z < num_cells; cell_z++) {
                int cell_index = cell_x + cell_y * num_cells + cell_z * num_cells * num_cells;

                for (const auto& offset : neighbor_offsets) {
                    int neighbor_x = (cell_x + std::get<0>(offset) + num_cells) % num_cells;
                    int neighbor_y = (cell_y + std::get<1>(offset) + num_cells) % num_cells;
                    int neighbor_z = (cell_z + std::get<2>(offset) + num_cells) % num_cells;
                    int neighbor_index = neighbor_x + neighbor_y * num_cells + neighbor_z * num_cells * num_cells;

                    // Add checks for neighbor_index
                    if (neighbor_index < 0 || neighbor_index >= total_cells) {
                        std::cerr << "ERROR: neighbor_index out of bounds: " << neighbor_index << std::endl;
                        continue; // Skip this neighbor
                    }

                    // Add check that cell_index is valid
                     if (cell_index < 0 || cell_index >= total_cells) {
                        std::cerr << "ERROR: cell_index out of bounds: " << cell_index << std::endl;
                        continue; // Skip this neighbor
                    }

                    for (int i : cell_list[cell_index]) {
                        for (int j : cell_list[neighbor_index]) {
                            if (i >= j) continue;

                            double dx = particles[j].x - particles[i].x;
                            double dy = particles[j].y - particles[i].y;
                            double dz = particles[j].z - particles[i].z;
                            dx -= L * round(dx / L);
                            dy -= L * round(dy / L);
                            dz -= L * round(dz / L);

                            double r = sqrt(dx * dx + dy * dy + dz * dz);

                            // Check for zero distance
                            if (r < 1e-12) {
                                std::cerr << "ERROR: Zero distance encountered between particles " << i << " and " << j << std::endl;
                                continue; // Skip this pair
                            }

                            if (r < rc) {
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
            }
        }
    }
}

// Function to update positions
void update_positions(int n, double L, double dt, Particle* particles, double* fx, double* fy, double* fz) {
    for (int i = 0; i < n; i++) {
        particles[i].x += particles[i].vx * dt + 0.5 * fx[i] * dt * dt / mass;
        particles[i].y += particles[i].vy * dt + 0.5 * fy[i] * dt * dt / mass;
        particles[i].z += particles[i].vz * dt + 0.5 * fz[i] * dt * dt / mass;
        // Apply periodic boundary conditions correctly
        particles[i].x = fmod(particles[i].x + L, L);
        particles[i].y = fmod(particles[i].y + L, L);
        particles[i].z = fmod(particles[i].z + L, L);

        //Clamp after position update
        particles[i].x = std::min(particles[i].x, L - 1e-12);
        particles[i].y = std::min(particles[i].y, L - 1e-12);
        particles[i].z = std::min(particles[i].z, L - 1e-12);
    }
}

// Function to update velocities
void update_velocities(int n, double dt, Particle* particles, double* fx, double* fy, double* fz, double* fx1, double* fy1, double* fz1) {
    for (int i = 0; i < n; i++) {
        particles[i].vx += 0.5 * (fx[i] + fx1[i]) * dt / mass;
        particles[i].vy += 0.5 * (fy[i] + fy1[i]) * dt / mass;
        particles[i].vz += 0.5 * (fz[i] + fz1[i]) * dt / mass;
    }
}

// Function to run the simulation
void run_simulation(int n, int steps, double dt, double L, double rc, int output_freq, Particle* particles, std::string basename) {
    double *fx = new double[n];
    double *fy = new double[n];
    double *fz = new double[n];
    double *fx1 = new double[n];
    double *fy1 = new double[n];
    double *fz1 = new double[n];

    double potential = 0.0, kinetic = 0.0, r = 0.0;

    std::cout << "Running simulation..." << std::endl;

    std::vector<std::vector<int>> cell_list(total_cells);

    initialize_positions(n, L, particles);
    initialize_velocities(n, particles);

    std::string energy_filename = basename + "_energy.dat";
    std::string trajectory_filename = basename + "_trajectory.xyz";

    std::ofstream energy_file(energy_filename);
    std::ofstream trajectory_file(trajectory_filename);

    
    

    build_cell_list(n, L, particles, cell_list);
    calculate_forces_linked_cell(n, L, particles, fx, fy, fz, cell_list);

    for (int step = 0; step < steps; step++) {

        for (int i = 0; i < n; i++) {
            fx1[i] = fx[i];
            fy1[i] = fy[i];
            fz1[i] = fz[i];
        }

        update_positions(n, L ,dt, particles, fx, fy, fz);
        build_cell_list(n, L ,particles, cell_list);
        calculate_forces_linked_cell(n, L, particles, fx, fy, fz, cell_list);
        update_velocities(n, dt, particles, fx, fy, fz, fx1, fy1, fz1);

        if (step % output_freq == 0) {
            potential = 0.0;
            for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    r = distance(particles[i].x, particles[i].y, particles[i].z, particles[j].x, particles[j].y, particles[j].z, L);
                    if (r < rc) {
                        potential += lj_potential(r);
                    }
                }
            }
            kinetic = 0.0;
            for (int i = 0; i < n; i++) {
                kinetic += 0.5 * mass * (particles[i].vx * particles[i].vx + particles[i].vy * particles[i].vy + particles[i].vz * particles[i].vz);
            }
            energy_file << step << " " << potential + kinetic << std::endl;
            
            trajectory_file << n << "\n\n";
            for (int i = 0; i < n; i++) {
                trajectory_file << "H " << particles[i].x << " " << particles[i].y << " " << particles[i].z << "\n";
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


int main() {
    Particle* particles = new Particle[N];
    
    std::string outdir = "run/";
    std::string basename = outdir + "linked_cell_3d";

    auto start = std::chrono::high_resolution_clock::now();
    run_simulation(N, steps, dt, L, rc, output_freq, particles, basename);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "Simulation Completed in: " << elapsed_seconds.count() << "s" << std::endl;

    delete[] particles;
    return 0;
}
