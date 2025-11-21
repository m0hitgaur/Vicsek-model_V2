import matplotlib.pyplot as plt
import numpy as np
import glob
from itertools import product
import os
import matplotlib.animation as animation

current_directory = os.getcwd()


flag=0
     


while(flag==0):

    check=input("Single frame(s/S) or animation(a/A)")

    if(check=="A" or check=="a"):

        # Parameters
        
        noise=input("Noise : ")
        angle=input("Angle : ")
        trial=input("Trial : ")
        

        detail=[]
        with open(current_directory+f'/data/config_data/parameters.txt', 'r') as file:
            for line in file:
                detail.append(line)
        Length_of_box = float(detail[1])   # Size of the grid
        numberoftrial=int(detail[10])
        maxiter=int(detail[3])
        Number_of_agents=int(detail[0])
        Lx=int(detail[1])
        Ly=float(detail[2])
        v0=float(detail[4])
        dt=float(detail[5])
        density=float(detail[8])
        rc=float(detail[9])


        times=[]
        for i in range( maxiter):
            if(i<10):tf=1
            if(i>10):tf=10
            if(i>100):tf=50
            if(i>1000):tf=100
            if(i%tf==0):times.append(i)


        # Set up the figure and axis
        fig, ax = plt.subplots()
        ax.set_aspect('equal')
        quiver = ax.quiver(0, 0, 0, 0, angles='xy', scale_units='xy', scale=1)



        # Function to load data from files
        def load_data(path):
            data = []
            with open(path, 'r') as f:
                    for line in f:
                        data.append(float(line.strip()))
            return np.array(data)
        colour=['r','b','g','m']




        def update(t):
            global quiver
            quiver.remove()  # Remove the previous quiver 

            path_template_x = current_directory+f'/Angle_{angle}/Noise_{noise}/flockingdata/positionx_{trial}_{times[t]}_.dat'
            path_template_y = current_directory+f'/Angle_{angle}/Noise_{noise}/flockingdata/positiony_{trial}_{times[t]}_.dat'
            path_template_theta =current_directory+f'/Angle_{angle}/Noise_{noise}/flockingdata/theta_{trial}_{times[t]}_.dat'
            
            # Load position and angle data
            positionx = load_data(path_template_x)
            positiony = load_data(path_template_y)
            theta = load_data(path_template_theta)
            
            detail=load_data(current_directory+f'/Angle_{angle}/Noise_{noise}/parameters.txt')
            grid_size = detail[1]   # Size of the grid
            
            # Extract data for the specific time step
            px = positionx
            py = positiony
            th = theta
            # Compute velocity components
            vx = np.cos(th)
            vy = np.sin(th)
            
            # Plot the quiver plot
            ax.set_xlim(Lx)
            ax.set_ylim(Ly)
            ax.set_title("N="+str(Number_of_agents)+r"|$\eta=$"+str(noise)+r"|$\alpha=$"+str(angle)+r"|$t=$"+str(times[t]))
            quiver=ax.quiver(px, py, vx, vy,angles='xy', scale_units='xy', scale=1)
            
            return quiver




            

        # Create the animation
        ani = animation.FuncAnimation(fig, update, frames=len(times), blit=False)

        # Show the animation
        plt.show()
        flag=1


    elif(check=="s" or check=="S"):
        # Parameters
        
        noise=input("Noise : ")
        angle=input("Angle : ")
        time=input("Time : ")
        trial=input("Trial : ")
        detail=[]
        with open(current_directory+f'/Angle_{angle}/Noise_{noise}/parameters.txt', 'r') as file:
            for line in file:
                detail.append(line)
        Length_of_box = float(detail[1])   # Size of the grid
        numberoftrial=int(detail[10])
        maxiter=int(detail[3])
        Number_of_agents=int(detail[0])
        Lx=int(detail[1])
        Ly=float(detail[2])
        v0=float(detail[4])
        dt=float(detail[5])
        density=float(detail[8])
        rc=float(detail[9])



        # Function to load data from files
        def load_data(path):
            data = []
            with open(path, 'r') as f:
                    for line in f:
                        data.append(float(line.strip()))
            return np.array(data)
        colour=['r','b','g','m']




        path_template_x = current_directory+f'/Angle_{angle}/Noise_{noise}/flockingdata/positionx_{trial}_{time}_.dat'
        path_template_y = current_directory+f'/Angle_{angle}/Noise_{noise}/flockingdata/positiony_{trial}_{time}_.dat'
        path_template_theta =current_directory+f'/Angle_{angle}/Noise_{noise}/flockingdata/theta_{trial}_{time}_.dat'
        
        # Load position and angle data
        positionx = load_data(path_template_x)
        positiony = load_data(path_template_y)
        theta = load_data(path_template_theta)
        
        detail=load_data(current_directory+f'/Angle_{angle}/Noise_{noise}/parameters.txt')
        grid_size = detail[1]   # Size of the grid
        
        # Extract data for the specific time step
        px = positionx
        py = positiony
        th = theta
        # Compute velocity components
        vx = np.cos(th)
        vy = np.sin(th)
        
        # Plot the quiver plot

        plt.quiver(px, py, vx, vy,angles='xy', scale_units='xy', scale=1)
        
        
        # Show the animation
        plt.show()
        flag=1

    else:
        print("Wrong option Try again!!")
        
