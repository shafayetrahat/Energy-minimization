#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <vector>

// Constants
double sigma = 1.0;
double epsilon = 1.0;
double mass = 1.0;
int N = 40; // number of particles
int steps = 1e6; // number of steps
double dt = 1e-4; // time step
double L = 10.0; // box length
double rc = 2.4; // cut-off distance
int output_freq = 200; // output frequency
char lattice_type = 't'; // lattice type ('t' for triangular, 's' for square)

// Linked Cell Algorithm parameters
double cell_size = rc + 0.1; // cell size equal to cutoff distance
int num_cells = static_cast<int>(L / cell_size); // number of cells in each dimension
int total_cells = num_cells * num_cells; // total number of cells

struct Particle {
    double x, y, vx, vy;
};

// Predefined neighboring cell offsets
const std::vector<std::pair<int, int>> neighbor_offsets = {
    {-1, -1}, {-1, 0}, {-1, 1},
    {0, -1}, {0, 0}, {0, 1},
    {1, -1}, {1, 0}, {1, 1}
};

// Function to initialize positions in a lattice
void initialize_positions(int n, double L, Particle* particles, char lattice_type) {
    if (lattice_type == 't') {
        int nx = ceil(sqrt(n));
        int ny = (n + nx - 1) / nx;
        double dx = L / nx;
        double dy = L / ny;
        int i = 0;
        for (int iy = 0; iy < ny; iy++) {
            for (int ix = 0; ix < nx; ix++) {
                if (i < n) {
                    particles[i].x = (ix + 0.5 * iy) * dx;
                    particles[i].y = iy * dy;
                    i++;
                }
            }
        }
    } else if (lattice_type == 's') {
        int nx = ceil(sqrt(n));
        int ny = (n + nx - 1) / nx;
        double dx = L / nx;
        double dy = L / ny;
        int i = 0;
        for (int iy = 0; iy < ny; iy++) {
            for (int ix = 0; ix < nx; ix++) {
                if (i < n) {
                    particles[i].x = ix * dx;
                    particles[i].y = iy * dy;
                    i++;
                }
            }
        }
    } else {
        std::cerr << "Invalid lattice type. Please enter 't' for triangular lattice or 's' for square lattice." << std::endl;
    }
}

// Function to initialize velocities
void initialize_velocities(int n, Particle* particles) {
    if (n <= 0) return;

    double vxcm = 0.0, vycm = 0.0;

    for (int i = 0; i < n; i++) {
        // Random velocities between -0.1 and 0.1
        particles[i].vx = 0.2 * (rand() / (RAND_MAX + 1.0) - 0.5);
        particles[i].vy = 0.2 * (rand() / (RAND_MAX + 1.0) - 0.5);
        

        vxcm += particles[i].vx;
        vycm += particles[i].vy;
    }

    vxcm /= n;
    vycm /= n;

    for (int i = 0; i < n; i++) {
        particles[i].vx -= vxcm;
        particles[i].vy -= vycm;
    }
}


// Function to calculate distance with PBC
double distance(double x1, double y1, double x2, double y2, double L) {
    double dx = x2 - x1;
    double dy = y2 - y1;
    dx -= L * round(dx / L);
    dy -= L * round(dy / L);
    return sqrt(dx * dx + dy * dy);
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
    
    for (auto& cell : cell_list) {
        cell.clear();  // Clear the cell list
    }

    int cell_x = 0, cell_y = 0, cell_index = 0;

    for (int i = 0; i < n; i++) {
        // Apply periodic boundary conditions correctly
        particles[i].x = fmod(particles[i].x + L, L);
        particles[i].y = fmod(particles[i].y + L, L);

        // Calculate cell indices
        cell_x = static_cast<int>(particles[i].x / cell_size);
        cell_y = static_cast<int>(particles[i].y / cell_size);

        // Ensure cell indices are within bounds
        if (cell_x >= num_cells) cell_x = num_cells - 1;
        if (cell_y >= num_cells) cell_y = num_cells - 1;

        cell_index = cell_x + cell_y * num_cells;
        cell_list[cell_index].push_back(i);
    }
}


// Function to calculate forces using the linked cell algorithm
void calculate_forces_linked_cell(int n, double L, double rc, Particle* particles, double* fx, double* fy, const std::vector<std::vector<int>>& cell_list) {
    std::fill(fx, fx + n, 0.0);
    std::fill(fy, fy + n, 0.0);


    int cell_x = 0, cell_y = 0, cell_index = 0, neighbor_x = 0, neighbor_y = 0, neighbor_index = 0;
    double r = 0.0, dx = 0.0, dy = 0.0, f = 0.0;

    for (cell_x = 0; cell_x < num_cells; cell_x++) {
        for (cell_y = 0; cell_y < num_cells; cell_y++) {
            cell_index = cell_x + cell_y * num_cells;

            for (const auto& offset : neighbor_offsets) {
                neighbor_x = (cell_x + offset.first + num_cells) % num_cells;
                neighbor_y = (cell_y + offset.second + num_cells) % num_cells;
                neighbor_index = neighbor_x + neighbor_y * num_cells;

                for (int i : cell_list[cell_index]) {
                    for (int j : cell_list[neighbor_index]) {
                        if (i >= j) continue;  // Avoid double-counting and self-interaction

                        r = distance(particles[i].x, particles[i].y, particles[j].x, particles[j].y, L);
                        if (r > 0 && r < rc) {
                            f = lj_force(r);
                            dx = particles[j].x - particles[i].x;
                            dy = particles[j].y - particles[i].y;
                            dx -= L * round(dx / L);
                            dy -= L * round(dy / L);

                            fx[i] += f * dx / r;
                            fy[i] += f * dy / r;
                            fx[j] -= f * dx / r;
                            fy[j] -= f * dy / r;
                        }
                    }
                }
            }
        }
    }
}


// Function to update positions
void update_positions(int n, double L, double dt, Particle* particles, double* fx, double* fy) {
    for (int i = 0; i < n; i++) {
        particles[i].x += particles[i].vx * dt + 0.5 * fx[i] * dt * dt / mass;
        particles[i].y += particles[i].vy * dt + 0.5 * fy[i] * dt * dt / mass;

        // Apply periodic boundary conditions correctly
        particles[i].x = fmod(particles[i].x + L, L);
        particles[i].y = fmod(particles[i].y + L, L);
    }
}


// Function to update velocities
void update_velocities(int n, double dt, Particle* particles, double *fx, double *fy, double *fx1, double *fy1) {
    for (int i = 0; i < n; i++) {
        particles[i].vx += 0.5 * (fx[i] + fx1[i]) * dt / mass;
        particles[i].vy += 0.5 * (fy[i] + fy1[i]) * dt / mass;
    }
}

// Function to run the simulation
void run_simulation(int n, int steps, double dt, double L, double rc, int output_freq, int seed, Particle* particles, char lattice_type, std::string basename) {
    srand(seed);
    double *fx = new double[n];
    double *fy = new double[n];
    double *fx1 = new double[n];
    double *fy1 = new double[n];

    double potential = 0.0, kinetic = 0.0, r = 0.0;

    std::cout << "Running simulation..."<< std::endl;

    std::vector<std::vector<int>> cell_list(total_cells);

    initialize_positions(n, L, particles, lattice_type);
    initialize_velocities(n, particles);

    std::string energy_filename = basename + "_energy.dat";
    std::string trajectory_filename = basename + "_trajectory.xyz";
    std::ofstream energy_file(energy_filename);
    std::ofstream trajectory_file(trajectory_filename);
   

    build_cell_list(n, L, particles, cell_list);
    calculate_forces_linked_cell(n, L, rc, particles, fx, fy, cell_list);

    for (int step = 0; step < steps; step++) {
        for (int i = 0; i < n; i++) {
            fx1[i] = fx[i];
            fy1[i] = fy[i];
        }

        update_positions(n, L, dt, particles, fx, fy);
        build_cell_list(n, L, particles, cell_list);
        calculate_forces_linked_cell(n, L, rc, particles, fx, fy, cell_list);
        update_velocities(n, dt, particles, fx, fy, fx1, fy1);

        if (step % output_freq == 0) {
            potential = 0.0;
            for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    r = distance(particles[i].x, particles[i].y, particles[j].x, particles[j].y, L);
                    if (r < rc) {
                        potential += lj_potential(r);
                    }
                }
            }
            kinetic = 0.0;
            for (int i = 0; i < n; i++) {
                kinetic += 0.5 * mass * (particles[i].vx * particles[i].vx + particles[i].vy * particles[i].vy);
            }
            energy_file << step << " " << potential + kinetic << std::endl;

            trajectory_file << n << "\n\n";
            for (int i = 0; i < n; i++) {
                trajectory_file << "H " << particles[i].x << " " << particles[i].y << " 0.0\n";
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

int main() {
    Particle* particles = new Particle[N];

    std::string outdir = "run/";
    std::string basename = outdir + "linked_cell_2d";

    // Start timing
    auto start = std::chrono::high_resolution_clock::now();

    run_simulation(N, steps, dt, L, rc, output_freq, 12345, particles, lattice_type, basename);

    // End timing
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    std::cout << "Simulation completed in " << elapsed.count() << " seconds." << std::endl;

    delete[] particles;

    return 0;
}