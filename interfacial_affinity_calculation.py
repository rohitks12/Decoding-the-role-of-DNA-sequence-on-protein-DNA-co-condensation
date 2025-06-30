#####################################################################################################################################################
#Here, we are identifying the interface of the condensates by using the probability profiles of the monomers being inside the condensate.
#The code is used in analysis for figure 5 (G) and figure 6 (C and D).
#Firstly, we identify the increasing trends (left interfaces) and decreasing trends (right interfaces) in the probability profile.
#Next, we compute the average interfacial affinity based on the probability of monomers in the defined interface and their interaction strength with the proteins.
#Finally, we plot the interfacial affinity along with probability and average interfacial affinity.
####################################################################################################################################################

import math
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


#affinity mapping as a string:
# for homogeneous DNA 2:2.00” means type 2 → 2.00 kT
# for heterogeneous DNA I 1:1.75,2:2.25” means type 1 → 1.75 kT, type 2 → 2.25 kT
#affinity_str = "2:2.00" #for homogeneous DNA
affinity_str = "1:1.75,2:2.25" #for heterogeneous DNA I


aff_map = {}
for pair in affinity_str.split(","):
    t, a = pair.split(":")
    aff_map[int(t)] = float(a)


aff_position = pd.read_csv(
    "/path/to/directory/TypeID_hetero_I.dat",
    sep=",", header=None,
    names=["particle_id", "types"]
)


aff_position["particle_id"] = range(len(aff_position))


aff_position["aff"] = aff_position["types"].map(aff_map)

print(aff_position.head())


directory = '/path/to/directory/'





def find_probability_windows(series, low_thresh=0.1, high_thresh=0.99):
    """
    We find asending windows in the probability series indicating the left interface
    first the window starts with probability 0.1 or less
    Next, code checks if the probability of i is more than average of the probabilities of the previous windows
    but less than 0.99.
    Once probability is equal to or more than 0.99, break.
    """
    windows = []
    n = len(series)
    i = 0

    while i < n:
        
        if series[i] <= low_thresh:
            start = i
            current = series[i]
            window = [current]
            prev_vals = [current]

            
            while i + 1 < n:
                avg_prev = np.mean(prev_vals)
                next_val = series[i + 1]

                if next_val > avg_prev:
                    i += 1
                    current = next_val
                    window.append(current)
                    prev_vals.append(current)
                else:
                    break

                
                if current >= high_thresh:
                    windows.append((start, i, window))
                    break
        i += 1

    return windows



def find_reverse_probability_windows(series, window_size=5, threshold=0.1):
    """
    We find the decending windows in the probability series indicating the right interface of the condensate.
    First, start with probability of more than or equal to 0.99.
    To capture the decreasing trend, see if the next value is less than the average of the values present in the window.
    Break when value is equal to or less than 0.1
    """
    windows = []  
    n = len(series)
    i = 0

    while i < n:
        if series[i] >= 0.99:  
            start = i  
            current_value = series[i]
            window = [current_value]
            prev_values = [current_value]

            while i + 1 < n:
                avg_prev = np.mean(prev_values)
                
                if series[i + 1] < avg_prev:  
                    prev_values.append(series[i + 1])
                    window.append(series[i + 1])
                    current_value = series[i + 1]
                else:
                    break

                if current_value <= threshold:
                    windows.append((start, i + 1, window))  
                    break

                i += 1  

        i += 1  

    return windows


# Looping through each probability file in the directory
for filename in os.listdir(directory):
    if filename.endswith(".dat"):  # filtering for .dat files only
        parts = filename.split('_')
        conc = parts[0]  
        replicate = parts[2].split('.')[0]  
        
        
        file_path = os.path.join(directory, filename)
        df = pd.read_csv(file_path, sep=",", header=0, names=['particle_id', 'probability'])
        probability = df['probability']

        
        windows = find_probability_windows(probability)
        print(windows)
        
        if windows:
            # Create transition data
            transition_data = []
            for start_idx, end_idx, prob_window in windows:
                particle_ids = list(range(start_idx, end_idx + 1))

                df_window = pd.DataFrame({
                    'particle_id': particle_ids,
                    'probability': prob_window
                })

                df_window = df_window.merge(aff_position, on='particle_id', how='left')
                transition_data.append(df_window)

           
            transition_df = pd.concat(transition_data, ignore_index=True)

           
            print(transition_df)

           
            transition_df['weighted_aff'] = transition_df['aff'] * transition_df['probability']
            transition_df['average_interfacial_affinity'] = transition_df['weighted_aff'].sum() / transition_df['probability'].sum()

            
            name_transition = f"Conc = {conc}, Replicate = {replicate}, Average Affinity = {transition_df['average_interfacial_affinity'].iloc[0]}"

            # Plotting the transition 0 to 1 (Affinity as bars and Probability as line)
            fig, ax1 = plt.subplots(figsize=(10, 6), dpi=300)
            ax1.bar(transition_df['particle_id'], transition_df['aff'], color='tab:blue', width=0.8)
            ax1.set_xlabel('Particle ID', fontsize=21)
            ax1.set_ylabel('Affinity', color='tab:blue', fontsize=21)
            ax1.tick_params(axis='y', labelcolor='tab:blue', labelsize=16)

            ax2 = ax1.twinx()
            ax2.plot(transition_df['particle_id'], transition_df['probability'], color='black', linewidth=2, label='Probability')
            ax2.set_ylabel('Probability', color='black', fontsize=21)
            ax2.tick_params(axis='y', labelcolor='black', labelsize=16)

            plt.title(f'Transition Plot (0 -> 1): Affinities (Bar) and Probabilities (Line)\n{name_transition}', fontsize=16)
            ax1.grid(True)

            
            plt.savefig(f'/path/to/directory/{conc}_{replicate}_transition_plot_0_to_1.svg', format='svg')
            plt.close()

        else:
            print("No transition region found in the probability array.")

        # Computing and plot the reverse transition (1 to 0)
        reverse_windows = find_reverse_probability_windows(probability)

        if reverse_windows:
            reverse_transition_data = []
            for start_idx, end_idx, prob_window in reverse_windows:
                particle_ids = list(range(start_idx, end_idx + 1))

                df_window = pd.DataFrame({
                    'particle_id': particle_ids,
                    'probability': prob_window
                })

                df_window = df_window.merge(aff_position, on='particle_id', how='left')
                reverse_transition_data.append(df_window)

            reverse_transition_df = pd.concat(reverse_transition_data, ignore_index=True)

            
            print(reverse_transition_df)

            
            reverse_transition_df['weighted_aff'] = reverse_transition_df['aff'] * reverse_transition_df['probability']
            reverse_transition_df['average_interfacial_affinity'] = reverse_transition_df['weighted_aff'].sum() / reverse_transition_df['probability'].sum()

            
            name_reverse = f"Conc = {conc}, Replicate = {replicate}, Average Affinity = {reverse_transition_df['average_interfacial_affinity'].iloc[0]}"

            # Plot the reverse transition 1 -> 0 (Affinity as bars and Probability as line)
            fig, ax1 = plt.subplots(figsize=(10, 6), dpi=300)
            ax1.bar(reverse_transition_df['particle_id'], reverse_transition_df['aff'], color='tab:blue', width=0.8)
            ax1.set_xlabel('Particle ID', fontsize=21)
            ax1.set_ylabel('Affinity', color='tab:blue', fontsize=21)
            ax1.tick_params(axis='y', labelcolor='tab:blue', labelsize=16)

            ax2 = ax1.twinx()
            ax2.plot(reverse_transition_df['particle_id'], reverse_transition_df['probability'], color='black', linewidth=2, label='Probability')
            ax2.set_ylabel('Probability', color='black', fontsize=21)
            ax2.tick_params(axis='y', labelcolor='black', labelsize=16)

            plt.title(f'Reverse Transition Plot (1 -> 0): Affinities (Bar) and Probabilities (Line)\n{name_reverse}', fontsize=16)
            ax1.grid(True)

            
            plt.savefig(f'/path/to/directory/{conc}_{replicate}_reverse_transition_plot_1_to_0.svg', format='svg')
            plt.close()

        else:
            print("No reverse transition region found in the probability array.")

