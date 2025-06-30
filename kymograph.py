################################################################################################################################################################################################
#The code plots a kymograph of the simulation based on distribution of DNA and proteins on the longitudinal axis
#To do this, the code first stores the position of all the particles on z-axis.
#Next, after trimming the distribution based on position of first and last monomer, it counts the number of particles in each bin
#where each bin is 2 sigma unit long.
#Plotting the kymograph and saving the data. 
################################################################################################################################################################################################

import os
import numpy as np
import pandas as pd
from natsort import natsorted
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

def process_configurations(directory, output_path):
    all_z_positions = []
    particle0_z = []  # For the first monomer position (left tethered end)
    particle499_z = []  # For the position of last monomer (right tethered end)

    # Processing each configuration file
    for filename in natsorted(os.listdir(directory)):
        f = os.path.join(directory, filename)
        print(f"Processing: {filename}")
        if os.path.isfile(f):
            try:
                
                df = pd.read_csv(f, sep=" ", header=None, 
                               names=["particle_id", "x_pos", "y_pos", "z_pos", "type"])
                df = df.drop([0, len(df)-1])  
                
                
                for col in df.columns:
                    df[col] = df[col].astype(str).str.replace(r'[\{\}\t]', '', regex=True)
                    df[col] = pd.to_numeric(df[col], errors='coerce')
                df = df.dropna()
                
                
                df = df.sort_values('particle_id')
                z_pos = df['z_pos'].values
                
                all_z_positions.append(z_pos)
                
                
                particle0_z.append(df[df['particle_id'] == 0]['z_pos'].values[0])
                particle499_z.append(df[df['particle_id'] == 499]['z_pos'].values[0])
                
            except Exception as e:
                print(f"Error processing {filename}: {str(e)}")
                continue

    if not all_z_positions:
        raise ValueError("No valid configuration files found in the directory")

    z_matrix = np.array(all_z_positions).T
    
    # computiong the z-axis range with padding of ±2 around the first and last monomer
    z_min = min(particle0_z) - 2  #first monomer position - 2
    z_max = max(particle499_z) + 2  # last monomer position + 2
    
    # setting bin width = 2 sigma
    bin_width = 2
    bins = np.arange(z_min, z_max + bin_width, bin_width)

    
    counts_matrix = np.array([np.histogram(z_matrix[:, t], bins=bins)[0] 
                           for t in range(z_matrix.shape[1])])
    normalized_counts = counts_matrix / z_matrix.shape[0]

    np.savez(output_path,
             normalized_counts=normalized_counts,
             bins=bins,
             xlim_min=z_min,
             xlim_max=z_max,
             y_ticks=[0, 2000, 4000, 6000, 8000, 10000],
             y_labels=['0', '1', '2', '3', '4', '5'])

def plot_kymograph(input_path, output_plot_path):
    """Generate kymograph plot from saved data"""
    data = np.load(input_path)
    normalized_counts = data['normalized_counts']
    bins = data['bins']
    xlim_min = data['xlim_min']
    xlim_max = data['xlim_max']

    fig, ax = plt.subplots(figsize=(10, 6), dpi=600)
    
    
    cax = ax.imshow(normalized_counts, 
                    cmap='gray', 
                    aspect='auto', 
                    origin='lower',
                    extent=[bins[0], bins[-1], 0, len(normalized_counts)],
                    interpolation='nearest')

   
    cbar = fig.colorbar(cax)
    cbar.set_ticks([0, 0.01, 0.02, 0.03])
    cbar.ax.tick_params(labelsize=25)
    cbar.set_label('Normalized frequency', fontsize=25)
    cax.set_clim(0.0, 0.03)

    
    ax.set_yticks(data['y_ticks'])
    ax.set_yticklabels(data['y_labels'], fontsize=25)
    
   
    ax.invert_yaxis()

    
    ax.set_xlabel('Z-position', fontsize=25)
    ax.set_ylabel(r'Time (x $10^7$ $\tau$)', fontsize=25)
    ax.set_xlim(xlim_min, xlim_max)
    ax.set_xticks([])

    
    fig.savefig(output_plot_path, format="svg", bbox_inches='tight', dpi=600)
    plt.close()  

if __name__ == "__main__":
    
    config_dir = '/path/to/output_files/'
    analysis_output = '/path/to/directory/kymo_data_0.6_rep_0.npz'
    process_configurations(config_dir, analysis_output)
    
    
    plot_output = '/path/to/directory/kymograph.svg'
    plot_kymograph(analysis_output, plot_output)
