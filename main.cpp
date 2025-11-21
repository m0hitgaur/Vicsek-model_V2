#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <algorithm>
#include <thread>
#include <filesystem>
using namespace std;
namespace fs = filesystem;

time_t trial_time,start_time=time(NULL) , finish_time;

bool create_directory(const string& path) {
    try {
        fs::create_directories(path);
        return fs::is_directory(path);
    } 
    catch (const fs::filesystem_error& e) {
        std::cerr << "Error creating directory '" << path << "': " << e.what() << '\n';
        return false;
    }
}
struct Particle {
    double x, y;
    double vx, vy;
    double ax, ay;
};

class Simulation {
private:
    vector<Particle> particles;
    double half_angle;
    int N;
    double Lx; 
    double Ly; 
    double dt;
    double v0;
    double noise;  // Noise strength
    double rc; // Cutoff radius
    double sigma,k;
    double alignment_strength;
    mt19937 gen;
    uniform_real_distribution<double> uniform_dist;
    
public:
    Simulation(int num_particles,double angle,double box_size_x,double box_size_y, double velo_mag,double timestep, double noise_strength, double align_str)
        : N(num_particles), half_angle(angle),Lx(box_size_x),Ly(box_size_y), dt(timestep), noise(noise_strength), alignment_strength(align_str),
           rc(3*sigma), k(1.0),sigma(1.0), uniform_dist(-1.0, 1.0),gen(12345) ,v0(velo_mag){
        particles.resize(N);
        initialize_particles();
    }
    
    void initialize_particles() {
        double spacing_x = Lx / sqrt(N);
        double spacing_y = Ly / sqrt(N);
        int idx = 0;
        
        for (int i = 0; i < sqrt(N) && idx < N; i++) {
            for (int j = 0; j < sqrt(N) && idx < N; j++) {
                
                particles[idx].x = i * spacing_x;
                particles[idx].y = j * spacing_y;
                
                double theta=M_PI*uniform_dist(gen);
                
                particles[idx].vx = v0*cos(theta);
                particles[idx].vy = v0*sin(theta);
                
                particles[idx].ax = 0;
                particles[idx].ay = 0;
                
                idx++;
            
            }
        }
    }       
    
    void pbc_posi(Particle & p){
        // Periodic boundary conditions
        if (p.x < 0) p.x += Lx;
        if (p.x >= Lx) p.x -= Lx;
        if (p.y < 0) p.y += Ly;
        if (p.y >= Ly) p.y -= Ly;
    }

    void pbc_velo(Particle & p){
        double theta=atan2(p.vy,p.vx);
        if(theta<-M_PI)theta= fmod(theta , M_PI)+M_PI;
        else if(theta>M_PI)theta= fmod(theta , M_PI)-M_PI;
        p.vx=v0 *cos(theta);
        p.vy=v0 *sin(theta);
    }
    double lennard_jones_force(double r) {
        const double rmin = 1e-3;          
        if (r > rc || r < rmin) return 0.0;
        double sr  = sigma / r;
        double sr6 = sr * sr * sr * sr * sr * sr;
        double sr12 = sr6 * sr6;
        return 24.0 * (2.0 * sr12 - sr6) / r;
    }
    double inter_particle_repulsive_force(double r) {
        const double rmin = 1e-3;
        if (r > 2*sigma || r < rmin) return 0.0;
        else return k*(r-2*sigma);
    }   
    double rij( Particle& p_i,  Particle& p_j) {
        double dx = p_i.x - p_j.x;
        double dy = p_i.y - p_j.y;

        if (dx > Lx/2) dx -= Lx;
        if (dx < -Lx/2) dx += Lx;
        if (dy > Ly/2) dy -= Ly;
        if (dy < -Ly/2) dy += Ly;
        double r=sqrt(dx*dx + dy*dy);
        if (r < 1e-3) r = 1e-3;
        return r;
    }
    
    void compute_forces() {
        for (Particle &p : particles) {
            p.ax = 0;
            p.ay = 0;
            }
        
        for (int i = 0; i < N; i++) {
            for (int j = i + 1; j < N; j++) {
                double dx = particles[j].x - particles[i].x;
                double dy = particles[j].y - particles[i].y;
                
                if (dx > Lx/2) dx -= Lx;
                if (dx < -Lx/2) dx += Lx;
                if (dy > Ly/2) dy -= Ly;
                if (dy < -Ly/2) dy += Ly;
                double r = sqrt(dx*dx + dy*dy );
                
                if (r < rc && r > 1e-6) {
                    double f = inter_particle_repulsive_force(r);
                    double fx = f * dx / r;
                    double fy = f * dy / r;
                    
                    particles[i].ax += fx;
                    particles[i].ay += fy;
                    
                    particles[j].ax -= fx;
                    particles[j].ay -= fy;
                    
                }
            }
        }
    }
    
    void velocity_alignment() {
    
        vector<double> avgx(N,0);
        vector<double> avgy(N,0),newtheta(N);
        vector<int> count(N,1); // starting from 1 as the particle itself is always counted
        
        for (int i=0;i<particles.size();i++) {
            double theta_i=atan2(particles[i].vy,particles[i].vx);
            avgx[i]+=cos(theta_i);   // just so that we can include the particle itself in average
            avgy[i]+=sin(theta_i);   // just so that we can include the particle itself in average
            
            for (int j = i+1; j < particles.size(); j++) 
            {   
                double theta_j=atan2(particles[j].vy,particles[j].vx);
                double dx=particles[j].x - particles[i].x;
                double dy=particles[j].y - particles[i].y; 
                if(dx>Lx/2) dx=dx-Lx;
                if(dx<-Lx/2) dx=dx+Lx ;  
                if(dy>Ly/2) dy=dy-Ly;
                if(dy<-Ly/2) dy=dy+Ly ;                 

                double rij = sqrt(pow(dx, 2) + pow(dy, 2));
                if (rij < 1e-4)rij=1e-4;
                double innerproduct_i=( (cos(theta_i) * (dx))+( sin(theta_i) * (dy) ) )/(rij); 
                double innerproduct_j= -1 * ( (cos(theta_j)*(dx) )+(sin(theta_j) * (dy)))/(rij); 
                        
                if(rij <= rc && innerproduct_i >= cos(half_angle) )      //////////////////////////// r_c or r_c/2??
                    {avgy[i]+= sin(theta_j);
                    avgx[i]+=cos(theta_j);
                    count[i]++;}
                
                if(rij <= rc && innerproduct_j >= cos(half_angle) ) 
                    {avgy[j]+= sin(theta_i);
                    avgx[j]+=cos(theta_i);
                    count[j]++;}
            }
            
        if(count[i]!=0)
        {avgx[i] /= static_cast<double>(count[i]);
        avgy[i]/=static_cast<double>(count[i]);}

        newtheta[i] = atan2(avgy[i],avgx[i]) + (uniform_dist(gen))*(noise/2);
        }
        
        for(int i=0;i<particles.size();i++){
            particles[i].vx += alignment_strength * v0* cos(newtheta[i]) * dt;
            particles[i].vy += alignment_strength * v0*sin(newtheta[i]) * dt;
            
        }    
    }

    void position_update(){
        for (auto& p : particles) {
            compute_forces();
            p.vx += p.ax;
            p.vy += p.ay ;
            pbc_velo(p);
            p.x += p.vx * dt;
            p.y += p.vy * dt ;           
            pbc_posi(p);
        }
    }

    void integrate() {
        position_update();
        velocity_alignment();
        
    }
    
    double velocity_order_parameter() {
        double sum_vx = 0, sum_vy = 0;
        
        for ( auto& p : particles) {
            sum_vx += p.vx;
            sum_vy += p.vy;    
        }
        
        double va = sqrt(sum_vx*sum_vx + sum_vy*sum_vy ) /(v0* N);
        return va;
    }
    
    void save_snapshot(int step,int trial) {
        ofstream file("data/config_data/config_" +to_string(trial)+"_" +to_string(step) + ".csv");
        for ( auto& p : particles) {
            file << p.x << "," << p.y << ","
                 << p.vx << "," << p.vy << "\n";
        }
        file.close();
    }
    
    void run_simulation(int tmax,int trialstart,int numberoftrials) {

        vector<int> times;
        for (int t = 0; t < tmax; t++) {
            bool should_record = false;
                if (t < 10) should_record = true;                    
                if (t<500 && t >= 100 && t % 10 == 0) should_record = true;   
                if (t >= 500&& t<1000 && t % 100 == 0) should_record = true;  
                if (t>=1000 && t % 500 == 0) should_record = true;                         
            if (should_record) times.push_back(t);
                }    
       
        vector<vector<double>> orderpara(numberoftrials-trialstart, vector <double>(times.size(),0) );
 
        for(int trial=trialstart;trial<numberoftrials;trial++){
                //ofstream f("/data/trial/")
                gen.seed(12345 + 10 * trial);
                string  path= "data/config_data/trial_"+ to_string(trial)+"/";
                create_directory( path);
                vector<double> order;
                trial_time=time(NULL);
                cout<<"\n"<<"Trial number : "<<trial<< " Out of "<<numberoftrials<<"    ";
            
        
            for (int t = 0; t < tmax; t++) {
                bool should_record = false;
                if (t < 10) should_record = true;                    
                if (t<500 && t >= 100 && t % 10 == 0) should_record = true;   
                if (t >= 500&& t<1000 && t % 100 == 0) should_record = true;  
                if (t>=1000 && t % 500 == 0) should_record = true;                        
                    
                
                integrate();
                
                if (should_record){
                    order.push_back(velocity_order_parameter());    
                    save_snapshot(t,trial);
                } 
                
                if (t % 100 == 0) cout  << t<<">>";
            } 
            
            orderpara[trial - trialstart]=order;  // check
            cout<<tmax<<"\n";    
            cout<<"Time to calculate trial = "  <<time(NULL)-trial_time<<" seconds  for   Angle : "+to_string(half_angle*180/M_PI)+" | Noise : "+to_string(noise)+" | Density : "+to_string(N/(Lx*Ly))+" | N = "+to_string(N)<<endl;  
            
            ofstream order_file("data/order_data/order_parameter_"+to_string(trial)+"_.csv");
            order_file<<"va,t\n";
            for (int i=0;i<order.size();i++)order_file<<order[i]<<","<<times[i]<<"\n";      
            order_file.close();
                    
            
        }    

        ofstream order_file("data/order_parameter.csv");
        string a="";
        for (int i=0;i<orderpara.size();i++)a+="trial_"+to_string(i)+",";      
        order_file<<a<<"t\n";
        for(int i=0;i<orderpara[0].size();i++){
            for(int j=0;j<orderpara.size();j++){
                order_file<<orderpara[j][i]<<",";
            }
            order_file<<times[i]<<"\n";
        }
        order_file.close();
        

        cout<<"\n"<<"Total time elapsed : "<< time(NULL) - start_time <<" seconds ";
        cout << "\nSimulation complete. Recorded " << times.size() << " snapshots." << endl;
    }
};

int main() {
    int N = 300;        // Number of particles
    double Lx = 20.0;    // Box size
    double Ly = 20.0;    // Box size
    double half_angle=M_PI;
    double v0=0.01;
    double dt = 0.01;   // Timestep
    double noise = 0.05; // Noise strength
    double align_str = 1;  // Alignment strength
    int tmax = 600;    // Maximum time
    int numberoftrials=1;
    int trialstart=0;
    Simulation sim(N, half_angle,Lx,Ly, v0,dt, noise, align_str);
    string path="data/order_data";
    create_directory(path);
    sim.run_simulation(tmax,trialstart,numberoftrials);
    
    return 0;
}