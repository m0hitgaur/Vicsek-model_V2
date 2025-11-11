#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <algorithm>

struct Particle {
    double x, y, z;
    double vx, vy, vz;
    double ax, ay, az;
};

class MDSimulation {
private:
    std::vector<Particle> particles;
    int N;
    double L;  // Box size
    double dt;
    double noise;  // Noise strength
    double rc; // Cutoff radius
    double sigma, epsilon;
    double alignment_strength;
    std::mt19937 rng;
    std::uniform_real_distribution<double> uniform_dist;
    
public:
    MDSimulation(int num_particles, double box_size, double timestep, double noise_strength, double align_str)
        : N(num_particles), L(box_size), dt(timestep), noise(noise_strength), 
          alignment_strength(align_str), rc(3.0), sigma(1.0), epsilon(1.0),
          rng(std::random_device{}()), uniform_dist(-1.0, 1.0) {
        
        particles.resize(N);
        initialize_particles();
    }
    
    void initialize_particles() {
        double spacing = L / std::cbrt(N);
        int idx = 0;
        
        for (int i = 0; i < std::cbrt(N) && idx < N; ++i) {
            for (int j = 0; j < std::cbrt(N) && idx < N; ++j) {
                for (int k = 0; k < std::cbrt(N) && idx < N; ++k) {
                    particles[idx].x = i * spacing;
                    particles[idx].y = j * spacing;
                    particles[idx].z = k * spacing;
                    
                    particles[idx].vx = uniform_dist(rng);
                    particles[idx].vy = uniform_dist(rng);
                    particles[idx].vz = uniform_dist(rng);
                    
                    particles[idx].ax = 0;
                    particles[idx].ay = 0;
                    particles[idx].az = 0;
                    
                    idx++;
                }
            }
        }
    }
    
    double lennard_jones_force(double r) {
        if (r > rc) return 0.0;
        double r6 = std::pow(r, 6);
        double r12 = r6 * r6;
        return 24.0 * epsilon * (2.0 * std::pow(sigma, 12) / r12 - std::pow(sigma, 6) / r6) / r;
    }
    
    double distance(const Particle& p1, const Particle& p2) {
        double dx = p1.x - p2.x;
        double dy = p1.y - p2.y;
        double dz = p1.z - p2.z;
        
        // Periodic boundary conditions
        if (dx > L/2) dx -= L;
        if (dx < -L/2) dx += L;
        if (dy > L/2) dy -= L;
        if (dy < -L/2) dy += L;
        if (dz > L/2) dz -= L;
        if (dz < -L/2) dz += L;
        
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    }
    
    void compute_forces() {
        for (auto& p : particles) {
            p.ax = 0;
            p.ay = 0;
            p.az = 0;
        }
        
        for (int i = 0; i < N; ++i) {
            for (int j = i + 1; j < N; ++j) {
                double dx = particles[j].x - particles[i].x;
                double dy = particles[j].y - particles[i].y;
                double dz = particles[j].z - particles[i].z;
                
                if (dx > L/2) dx -= L;
                if (dx < -L/2) dx += L;
                if (dy > L/2) dy -= L;
                if (dy < -L/2) dy += L;
                if (dz > L/2) dz -= L;
                if (dz < -L/2) dz += L;
                
                double r = std::sqrt(dx*dx + dy*dy + dz*dz);
                
                if (r < rc && r > 1e-6) {
                    double f = lennard_jones_force(r);
                    double fx = f * dx / r;
                    double fy = f * dy / r;
                    double fz = f * dz / r;
                    
                    particles[i].ax += fx;
                    particles[i].ay += fy;
                    particles[i].az += fz;
                    
                    particles[j].ax -= fx;
                    particles[j].ay -= fy;
                    particles[j].az -= fz;
                }
            }
        }
    }
    
    void velocity_alignment() {
        double avg_vx = 0, avg_vy = 0, avg_vz = 0;
        
        for (const auto& p : particles) {
            avg_vx += p.vx;
            avg_vy += p.vy;
            avg_vz += p.vz;
        }
        
        avg_vx /= N;
        avg_vy /= N;
        avg_vz /= N;
        
        for (auto& p : particles) {
            p.vx += alignment_strength * avg_vx * dt;
            p.vy += alignment_strength * avg_vy * dt;
            p.vz += alignment_strength * avg_vz * dt;
        }
    }
    
    void add_noise() {
        for (auto& p : particles) {
            p.vx += noise * uniform_dist(rng) * dt;
            p.vy += noise * uniform_dist(rng) * dt;
            p.vz += noise * uniform_dist(rng) * dt;
        }
    }
    
    void integrate() {
        compute_forces();
        
        for (auto& p : particles) {
            p.vx += p.ax * dt;
            p.vy += p.ay * dt;
            p.vz += p.az * dt;
            
            p.x += p.vx * dt;
            p.y += p.vy * dt;
            p.z += p.vz * dt;
            
            // Periodic boundary conditions
            if (p.x < 0) p.x += L;
            if (p.x >= L) p.x -= L;
            if (p.y < 0) p.y += L;
            if (p.y >= L) p.y -= L;
            if (p.z < 0) p.z += L;
            if (p.z >= L) p.z -= L;
        }
        
        velocity_alignment();
        add_noise();
    }
    
    double velocity_order_parameter() {
        double sum_vx = 0, sum_vy = 0, sum_vz = 0;
        
        for (const auto& p : particles) {
            sum_vx += p.vx;
            sum_vy += p.vy;
            sum_vz += p.vz;
        }
        
        double magnitude = std::sqrt(sum_vx*sum_vx + sum_vy*sum_vy + sum_vz*sum_vz) / N;
        return magnitude;
    }
    
    void save_snapshot(int step) {
        std::ofstream file("snapshots_" + std::to_string(step) + ".dat");
        for (const auto& p : particles) {
            file << p.x << " " << p.y << " " << p.z << " "
                 << p.vx << " " << p.vy << " " << p.vz << "\n";
        }
        file.close();
    }
    
    double velocity_correlation(int distance_bin, double bin_width) {
        double correlation = 0;
        int count = 0;
        
        for (int i = 0; i < N; ++i) {
            for (int j = i + 1; j < N; ++j) {
                double r = distance(particles[i], particles[j]);
                
                if (r >= distance_bin * bin_width && r < (distance_bin + 1) * bin_width) {
                    double v_dot = particles[i].vx * particles[j].vx +
                                  particles[i].vy * particles[j].vy +
                                  particles[i].vz * particles[j].vz;
                    
                    double v_mag_i = std::sqrt(particles[i].vx*particles[i].vx + 
                                              particles[i].vy*particles[i].vy + 
                                              particles[i].vz*particles[i].vz);
                    double v_mag_j = std::sqrt(particles[j].vx*particles[j].vx + 
                                              particles[j].vy*particles[j].vy + 
                                              particles[j].vz*particles[j].vz);
                    
                    if (v_mag_i > 1e-6 && v_mag_j > 1e-6) {
                        correlation += v_dot / (v_mag_i * v_mag_j);
                        count++;
                    }
                }
            }
        }
        
        return count > 0 ? correlation / count : 0.0;
    }
    
    void run_simulation(int tmax) {
        std::vector<int> times;
        std::ofstream order_file("order_parameter.dat");
        
        for (int t = 0; t < tmax; ++t) {
            int ti = 0;
            if (t < 10) ti = 0;
            else if (t < 100) ti = 10;
            else if (t < 500) ti = 100;
            else if (t < 1000) ti = 500;
            else ti = 1000;
            
            if (ti > 0 && t % ti == 0) {
                times.push_back(t);
            }
            
            integrate();
            
            double phi = velocity_order_parameter();
            order_file << t << " " << phi << "\n";
            
            if (t % 100 == 0) {
                std::cout << "Step " << t << " Order Parameter: " << phi << std::endl;
                save_snapshot(t);
            }
        }
        
        order_file.close();
        
        // Calculate correlation functions
        std::ofstream vel_corr_file("velocity_correlation.dat");
        std::ofstream conn_corr_file("connected_correlation.dat");
        
        double bin_width = 0.5;
        int max_bins = static_cast<int>(L / (2 * bin_width));
        
        for (int b = 0; b < max_bins; ++b) {
            double v_corr = velocity_correlation(b, bin_width);
            double r_center = (b + 0.5) * bin_width;
            vel_corr_file << r_center << " " << v_corr << "\n";
            
            // Connected correlation (simplified: same as velocity correlation here)
            conn_corr_file << r_center << " " << v_corr << "\n";
        }
        
        vel_corr_file.close();
        conn_corr_file.close();
        
        std::cout << "\nSimulation complete. Recorded " << times.size() << " snapshots." << std::endl;
    }
};

int main() {
    int N = 100;        // Number of particles
    double L = 20.0;    // Box size
    double dt = 0.01;   // Timestep
    double noise = 0.5; // Noise strength
    double align_str = 0.1;  // Alignment strength
    int tmax = 2000;    // Maximum time
    
    MDSimulation sim(N, L, dt, noise, align_str);
    sim.run_simulation(tmax);
    
    return 0;
}