import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import matplotlib.animation as animation

def get_recording_times(maxiter):
    """Generate list of time steps that were recorded"""
    times = []
    for i in range(maxiter):
        if(i<10):times.append(i)
        if( 100<=i and i<500 and i % 10 == 0):times.append(i)
        if(i>1000 and i % 100 == 0):times.append(i)

    return times

def load_simulation_parameters():
    """Load parameters from CSV file"""
    param_path = os.path.join(os.getcwd(), 'data', 'parameters.csv')
    df = pd.read_csv(param_path)
    return {
        'N': int(df["N"].values),
        'Lx': int(df["Lx"].values),
        'Ly': int(df["Ly"].values),
        'alpha': float(df["alpha"].values),
        'noise': float(df["eta"].values),
        'maxiter': int(df["maxiter"].values)
    }
size=2
def plot_single_frame():
    """Plot a single time frame"""
    time_step = input("Time: ")
    trial = input("Trial: ")
    
    params = load_simulation_parameters()
    
    # Load data
    data_path = os.path.join(os.getcwd(), 'data', 'config_data',f'trial_{trial}', f'config_{time_step}.csv')
    df = pd.read_csv(data_path)
    
    # Extract arrays (not single values!)
    px = df["x"].values
    py = df["y"].values
    vx = df["vx"].values
    vy = df["vy"].values
    
    # Create plot
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_aspect('equal')
    ax.set_xlim(0, params['Lx'])
    ax.set_ylim(0, params['Ly'])·
    
    ax.scatter(px, py,s=size)
    ax.set_title(f"N={params['N']} | η={params['noise']} | α={params['alpha']}° | t={time_step}")
    
    plt.show()

def plot_animation():
    """Create animation over time"""
    trial = input("Trial: ")
    
    params = load_simulation_parameters()
    
    times = get_recording_times(params['maxiter'])
    
    # Set up figure
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_aspect('equal')
    ax.set_xlim(0, params['Lx'])
    ax.set_ylim(0, params['Ly'])
    
    # Initialize empty quiver
    quiver = ax.scatter([], [])
    
    def update(frame_idx):
        nonlocal quiver
        t = times[frame_idx]
        # Load data
        data_path = os.path.join(os.getcwd(), 'data', 'config_data',f'trial_{trial}', f'config_{t}.csv')
        
        df = pd.read_csv(data_path)
        
        # Extract arrays
        px = df["x"].values
        py = df["y"].values
        vx = df["vx"].values
        vy = df["vy"].values
        
        # Update quiver
        quiver.remove()
        quiver = ax.scatter(px, py,s=size)
        ax.set_title(f"N={params['N']} | η={params['noise']} | α={params['alpha']}° | t={t}")
        
        return quiver,
    
    # Create animation
    ani = animation.FuncAnimation(fig, update, frames=len(times), blit=False, interval=50)
    plt.show()

# Main loop
while True:
    check = input("Single frame(s/S) or animation(a/A): ").strip().lower()
    
    if check == 'a':
        plot_animation()
        break
    elif check == 's':
        plot_single_frame()
        break
    else:
        print("Wrong option. Try again!!")