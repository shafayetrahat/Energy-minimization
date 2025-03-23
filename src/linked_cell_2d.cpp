#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <vector>

// Constants
constexpr double sigma = 1.0;
constexpr double epsilon = 1.0;
constexpr double mass = 1.0;
constexpr int N = 500; // number of particles
constexpr int steps = 1e6; // number of steps
constexpr double dt = 1e-4; // time step
constexpr double L = 50.0; // box length
constexpr double rc = 2.5; // cut-off distance
constexpr int output_freq = 100; // output frequency

// Linked Cell Algorithm parameters
constexpr double cell_size = rc; // cell size equal to cutoff distance
constexpr int num_cells = static_cast<int>(L / cell_size); // number of cells in each dimension
constexpr int total_cells = num_cells * num_cells; // total number of cells

// Function to initialize positions in a lattice
void initialize_positions(int n, double L, double *x, double *y, char lattice_type) {
    if (lattice_type == 't') {
        int nx = ceil(sqrt(n));
        int ny = (n + nx - 1) / nx;
        double dx = L / nx;
        double dy = L / ny;
        int i = 0;
        for (int iy = 0; iy < ny; iy++) {
            for (int ix = 0; ix < nx; ix++) {
                if (i < n) {
                    x[i] = (ix + 0.5 * iy) * dx;
                    y[i] = iy * dy;
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
                    x[i] = ix * dx;
                    y[i] = iy * dy;
                    i++;
                }
            }
        }
    } else {
        std::cerr << "Invalid lattice type. Please enter 't' for triangular lattice or 's' for square lattice." << std::endl;
    }
}

// Function to initialize velocities
void initialize_velocities(int n, double *vx, double *vy) {
    double vxcm = 0.0;
    double vycm = 0.0;
    for (int i = 0; i < n; i++) {
        vx[i] = 1e-1 * (2.0 * rand() / RAND_MAX - 1.0);
        vy[i] = 1e-1 * (2.0 * rand() / RAND_MAX - 1.0);
        vxcm += vx[i];
        vycm += vy[i];
    }
    vxcm /= n;
    vycm /= n;
    for (int i = 0; i < n; i++) {
        vx[i] -= vxcm;
        vy[i] -= vycm;
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
void build_cell_list(int n, double L, double *x, double *y, std::vector<int> *cell_list) {
    for (int i = 0; i < total_cells; i++) {
        cell_list[i].clear(); // Clear the cell list
    }

    for (int i = 0; i < n; i++) {
        // Ensure positions are within the box [0, L)
        x[i] = x[i] - L * floor(x[i] / L);
        y[i] = y[i] - L * floor(y[i] / L);

        // Calculate cell indices
        int cell_x = static_cast<int>(x[i] / cell_size);
        int cell_y = static_cast<int>(y[i] / cell_size);

        // Ensure cell indices are within bounds using modulo
        cell_x = (cell_x + num_cells) % num_cells;
        cell_y = (cell_y + num_cells) % num_cells;

        int cell_index = cell_x + cell_y * num_cells;

        if (cell_index >= 0 && cell_index < total_cells) {
            cell_list[cell_index].push_back(i);
        } else {
            std::cerr << "Error: Invalid cell index " << cell_index << " for particle " << i << std::endl;
        }
    }
}

// Function to calculate forces using the linked cell algorithm
void calculate_forces_linked_cell(int n, double L, double rc, double *x, double *y, double *fx, double *fy, std::vector<int> *cell_list) {
    for (int i = 0; i < n; i++) {
        fx[i] = 0.0;
        fy[i] = 0.0;
    }

    for (int cell_x = 0; cell_x < num_cells; cell_x++) {
        for (int cell_y = 0; cell_y < num_cells; cell_y++) {
            int cell_index = cell_x + cell_y * num_cells;

            // Loop over neighboring cells (including itself)
            for (int neighbor_x = -1; neighbor_x <= 1; neighbor_x++) {
                for (int neighbor_y = -1; neighbor_y <= 1; neighbor_y++) {
                    int neighbor_cell_x = (cell_x + neighbor_x + num_cells) % num_cells;
                    int neighbor_cell_y = (cell_y + neighbor_y + num_cells) % num_cells;
                    int neighbor_cell_index = neighbor_cell_x + neighbor_cell_y * num_cells;

                    // Loop over particles in the current cell
                    for (int i : cell_list[cell_index]) {
                        // Loop over particles in the neighboring cell
                        for (int j : cell_list[neighbor_cell_index]) {
                            if (i < j) { // Avoid double-counting
                                double r = distance(x[i], y[i], x[j], y[j], L);
                                if (r > 0 && r < rc) { // Ensure r is not zero and within cutoff
                                    double f = lj_force(r);
                                    double dx = x[j] - x[i];
                                    double dy = y[j] - y[i];
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
    }
}

// Function to update positions
void update_positions(int n, double L, double dt, double *x, double *y, double *vx, double *vy, double *fx, double *fy) {
    for (int i = 0; i < n; i++) {
        x[i] += vx[i] * dt + 0.5 * fx[i] * dt * dt / mass;
        y[i] += vy[i] * dt + 0.5 * fy[i] * dt * dt / mass;
        x[i] = x[i] - L * floor(x[i] / L); // Apply PBC
        y[i] = y[i] - L * floor(y[i] / L); // Apply PBC
    }
}

// Function to update velocities
void update_velocities(int n, double dt, double *vx, double *vy, double *fx, double *fy, double *fx1, double *fy1) {
    for (int i = 0; i < n; i++) {
        vx[i] += 0.5 * (fx[i] + fx1[i]) * dt / mass;
        vy[i] += 0.5 * (fy[i] + fy1[i]) * dt / mass;
    }
}

// Function to run the simulation
void run_simulation(int n, int steps, double dt, double L, double rc, int output_freq, int seed, double *x, double *y, double *vx, double *vy, char lattice_type, std::string basename) {
    srand(seed);
    double *fx = new double[n];
    double *fy = new double[n];
    double *fx1 = new double[n];
    double *fy1 = new double[n];
    std::vector<int> *cell_list = new std::vector<int>[total_cells];

    initialize_positions(n, L, x, y, lattice_type);
    initialize_velocities(n, vx, vy);

    std::string energy_filename = basename + "_energy.dat";
    std::string trajectory_filename = basename + "_trajectory.xyz";
    std::ofstream energy_file(energy_filename);
    std::ofstream trajectory_file(trajectory_filename);

    for (int step = 0; step < steps; step++) {
        build_cell_list(n, L, x, y, cell_list);
        calculate_forces_linked_cell(n, L, rc, x, y, fx, fy, cell_list);

        for (int i = 0; i < n; i++) {
            fx1[i] = fx[i];
            fy1[i] = fy[i];
        }

        update_positions(n, L, dt, x, y, vx, vy, fx, fy);
        calculate_forces_linked_cell(n, L, rc, x, y, fx, fy, cell_list);
        update_velocities(n, dt, vx, vy, fx, fy, fx1, fy1);

        if (step % output_freq == 0) {
            double energy = 0.0;
            for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    double r = distance(x[i], y[i], x[j], y[j], L);
                    if (r < rc) {
                        energy += lj_potential(r);
                    }
                }
            }
            double kinetic = 0.0;
            for (int i = 0; i < n; i++) {
                kinetic += 0.5 * mass * (vx[i] * vx[i] + vy[i] * vy[i]);
            }
            energy_file << step << " " << energy + kinetic << std::endl;

            trajectory_file << n << "\n\n";
            for (int i = 0; i < n; i++) {
                trajectory_file << "H " << x[i] << " " << y[i] << " 0.0\n";
            }
        }
    }

    energy_file.close();
    trajectory_file.close();
    delete[] fx;
    delete[] fy;
    delete[] fx1;
    delete[] fy1;
    delete[] cell_list;
}

int main() {
    double *x = new double[N];
    double *y = new double[N];
    double *vx = new double[N];
    double *vy = new double[N];
    char lattice_type = 't';
    std::string outdir = "run/";
    std::string basename = outdir + "link_2d";

    // Start timing
    auto start = std::chrono::high_resolution_clock::now();

    run_simulation(N, steps, dt, L, rc, output_freq, 12345, x, y, vx, vy, lattice_type, basename);

    // End timing
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    std::cout << "Simulation completed in " << elapsed.count() << " seconds." << std::endl;

    delete[] x;
    delete[] y;
    delete[] vx;
    delete[] vy;

    return 0;
}