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
    double x_new,y_new;
    double vx_new,vy_new;
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
    int trial;
    vector<int> times;
    vector<double> order;
    mt19937 gen;
    uniform_real_distribution<double> uniform_dist,uniform_dist_N;
    string folder_path;
public:
    Simulation(int num_particles,double angle,double box_size_x,double box_size_y,int trial_, double velo_mag,double timestep, double noise_strength,string folderpath)
        : N(num_particles), half_angle(angle),Lx(box_size_x),Ly(box_size_y), dt(timestep), noise(noise_strength),
           rc(3*sigma), k(1.0),sigma(1.0), uniform_dist(-1.0, 1.0),uniform_dist_N(0, N),gen(12345),v0(velo_mag),folder_path(folderpath),trial(trial_){
        particles.resize(N);
        initialize_particles();
        
    }
    void save_order_data(){
        ofstream order_file(folder_path+"order_data/order_parameter_"+to_string(trial)+"_.csv");
        order_file<<"va,t\n";
        for (int i=0;i<order.size();i++)order_file<<order[i]<<","<<times[i]<<"\n";      
        order_file.close();

    }

    vector<double> get_order_data(){return order;}
    vector<int> get_time_data(){return times;}

    void initialize_particles() {
        gen.seed(12345 + 10 * trial);
        double rho = N/(Lx*Ly);
        if(rho>=1){cout<<"Particles Intialize Failure ";}
        int id=0;    
        for(int j=0;j<Lx &&id<N;j++){
            for(int k=0;k<Ly&&id<N;k++){

                particles[id].x = j ;
                particles[id].y = k;
                particles[id].x_new = j ;
                particles[id].y_new = k;

                double theta=M_PI*uniform_dist(gen);
    
                particles[id].vx = v0*cos(theta);
                particles[id].vy = v0*sin(theta);
                particles[id].vx_new = v0*cos(theta);
                particles[id].vy_new = v0*sin(theta);
                particles[id].ax = 0;
                particles[id].ay = 0;
                id++;                
            } 
        }

            
        
    
        //for(int t=0;t<100;t++)position_update();
    }          
    double dot_product(double theta,double dx,double dy,double rij){
        return ( (cos(theta) * (dx))+( sin(theta) * (dy) ) )/(rij);
    }
    void pbc_posi(Particle & p){
        // Periodic boundary conditions
        if (p.x < 0) p.x = fmod(p.x,Lx) +Lx;
        if (p.x > Lx) p.x = fmod(p.x,Lx);
        if (p.y < 0) p.y = fmod(p.y,Ly) +Ly;
        if (p.y > Ly) p.y = fmod(p.y,Ly);
    }
    double minimum_image(double dx,double Lx){
        if(dx>Lx/2) dx=dx-Lx;
        if(dx<-Lx/2) dx=dx+Lx ;  
        return dx;
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
        dx=minimum_image(dx,Lx);
        dy=minimum_image(dx,Ly);
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
                dx=minimum_image(dx,Lx);
                dy=minimum_image(dy,Ly);
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
                double dx=minimum_image(particles[j].x - particles[i].x,Lx);
                double dy=minimum_image(particles[j].y - particles[i].y,Ly);                 
                double rij = sqrt(pow(dx, 2) + pow(dy, 2));
                double innerproduct_i=  dot_product(theta_i,dx,dy,rij); 
                double innerproduct_j= -1 * dot_product(theta_j,dx,dy,rij); 
                        
                if(rij <= rc && innerproduct_i >= cos(half_angle) )      
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
            particles[i].vx = v0* cos(newtheta[i]) ;
            particles[i].vy = v0*sin(newtheta[i]) ;
            
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
        velocity_alignment(); 
        position_update(); 
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
        ofstream file(folder_path+"config_data/trial_"+to_string(trial)+"/config_" +to_string(step) + ".csv");
        file<<"x,y,vx,vy\n";
        for ( auto& p : particles) {
            file << p.x << "," << p.y << ","
                 << p.vx << "," << p.vy << "\n";
        }
        file.close();
    }   

    void run_simulation(int tmax,int trialstart,int numberoftrials,vector<bool>time_record) {
        ofstream f(folder_path+"parameters.csv");
        string head="N,Lx,Ly,alpha,v0,dt,eta,maxiter,numberoftrials\n";
        f<< head;
        f<<N<<","<<Lx<<","<<Ly<<","<<half_angle<<","<<v0<<","<<dt<<","<<noise<<","<<tmax<<","<<numberoftrials; 
        f.close();
        int time_counter=0;
        for (int t = 0; t < tmax; t++) {
            integrate();
            if (time_record[t]){
                order.push_back(velocity_order_parameter());    
                save_snapshot(t,trial);
                times.push_back(t);
                } 
            if (t % 100 == 0) cout  << t<<">>";
        } 
        cout<<tmax<<"\n"; 
        save_order_data();
        cout << "\nSimulation complete. Recorded " << times.size() << " snapshots." << endl;
    }

};



void save_order(vector<vector<double>>orderpara,vector<int>times,int trialstart,string folder_path){        
    ofstream order_file(folder_path+"order_parameter.csv");
    string a="";
    for (int i=0;i<orderpara.size();i++)a+="trial_"+to_string(i+trialstart)+",";      
    order_file<<a<<"t\n";
    for(int i=0;i<orderpara[0].size();i++){
        for(int j=0;j<orderpara.size();j++){
            order_file<<orderpara[j][i]<<",";
        }
        order_file<<times[i]<<"\n";
    }
    order_file.close();
    }

int main() {
    int N = 100;        // Number of particles
    double Lx = 15.0;    // Box size
    double Ly = 15.0;    // Box size
    double half_angle=M_PI;
    double v0=0.01;
    double dt = 0.01;   // Timestep
    double noise = 0.05; // Noise strength
    int tmax = 2000;    // Maximum time
    int numberoftrials=1;
    int trialstart=0;
    time_t trial_time,start_time=time(NULL) , finish_time;
    vector<bool> time_record;    
    for (int t = 0; t < tmax; t++) {
        bool should_record = false;
            if (t <= 10) should_record = true;
            if(t>=10 && t<100&& t%10==0) should_record = true;                   
            if (t<1000 && t >= 100 && t % 50 == 0) should_record = true;   
            if (t>=1000 && t % 100 == 0) should_record = true;                        
            time_record.push_back(should_record);
            }
    
    vector <vector<double>> order_data(numberoftrials-trialstart);
    vector<int> times;
    string folder_path;
    for(int trial=trialstart;trial<numberoftrials;trial++){
        folder_path="data/";
        create_directory(folder_path+"config_data/trial_"+ to_string(trial)+"/");
        trial_time=time(NULL);
        cout<<"\n"<<"Trial number : "<<trial<< " Out of "<<numberoftrials<< " | Angle : "+to_string(half_angle*180/M_PI)+" | Noise : "+to_string(noise)+" | Density : "+to_string(N/(Lx*Ly))+" | N = "+to_string(N)<<" | "<<endl ;
        
        Simulation sim(N, half_angle,Lx,Ly,trial,v0,dt, noise,folder_path);
        sim.run_simulation(tmax,trialstart,numberoftrials,time_record);
        
        cout<<"Time to calculate trial = "  <<time(NULL)-trial_time<<" seconds ";  
        
        order_data[trial]=sim.get_order_data();       
        if(trial==numberoftrials-1)times=sim.get_time_data();
    }   
    
    save_order(order_data,times,trialstart,folder_path); 
    
    cout<<"\n"<<"Total time elapsed : "<< time(NULL) - start_time <<" seconds ";
    return 0;
}