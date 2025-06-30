#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 17:48:41 2024

@author: rohitks
"""
###############################################################################################################################################################
#The code calculates the force exerted on the bare DNA (DNA outside the condensate) due to protein-DNA co-condensation.
#This code has been used in analysis for figure 3 (A and B), figure 6 (A and B).
#Firstly, the code uses DBSCAN algorithm to identify the monomers which are part of the condensate.
#Next, it calculates the bond lengths of monomer pairs outside the condensate.
#Next, it computes the force using l_0 which is mean bond length of the DNA in absence of proteins.
#Finally, it converts force measured in kBT/sigma into pN.
###############################################################################################################################################################

#imporing libraries
import math
import os 
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.cluster import DBSCAN #clustering algorithm
from natsort import natsorted #for sorting the files in directory
from sklearn.neighbors import NearestNeighbors #for finding nearest neighbors 
from scipy.signal import savgol_filter
from scipy.signal import find_peaks

#for DBSCAN, given parameter Min_points = 2 * dimension [sander et al. 1998]
#calculate distance with 6th nearest neighbor given Min_points = 6
def calculate_nearest_neighbors(i):
    neighb = NearestNeighbors(n_neighbors=6)
    nbrs = neighb.fit(i)
    distances, indices = nbrs.kneighbors(i)
    return distances, indices

#arrange the distance between pairs in asending order
#as rcut for attractive LJ-potential is 2.5 sigma, any particle pair which has distance more than rcut is not interacting so it is removed from calculating eps parameter
def sort_and_filter_distances(distances, cutoff_value=2.5):
    distances = np.sort(distances, axis=0)
    distances = distances[:, 5]
    distances = distances[distances <= cutoff_value]
    return distances


#as the k-distance curve has noise, point of maximum curvature can go to unreasonable values on the k-distance curve
#here as a smoothing filter we are using savgol filter which fits a polynomial of given order to a given window of data points on the curve
def smooth_data(distances):
    x = range(len(distances))
    smoothed_data = savgol_filter(distances, window_length=199, polyorder=3, mode='nearest')
    return smoothed_data


#the distance point on the curve for which the second derivative (curvature) is maximum (knee point) is treated as epsilon 
def calculate_derivatives(smoothed_data):
    x = range(len(smoothed_data))
    slope = np.gradient(smoothed_data, x)
    second_derivative = np.gradient(slope, x)
    return slope, second_derivative

def find_curvature_peaks(second_derivative):
    peaks, _ = find_peaks(second_derivative)
    return peaks

def find_max_curvature_point(distances, second_derivative):
    max_curvature_index = np.argmax(second_derivative)
    max_curvature_point = distances[max_curvature_index]
    return max_curvature_index, max_curvature_point



#lists
bond_length = [] #list for bond length for each configurations
mean_bond_length = [] #list of mean bond lengths for all configurations 
epsilon_list = []




#volume of condensates
vol = []

# assigning directory
directory = '/home/rohitks/co-condensation/github/output_files'

for filename in natsorted(os.listdir(directory)):
	f = os.path.join(directory, filename)

	if os.path.isfile(f):
          i = pd.read_csv(f, sep=" ", header=None, names=["index","x_pos", "y_pos", "z_pos", "type"])
          print(f)
          
          i = i.drop(0, axis=0)
          b = len(i) 
          i = i.drop(b, axis=0)
          i = i.drop('type', axis=1)
          i = i.drop('index', axis=1)
          
          
          distances, indices = calculate_nearest_neighbors(i)
          filtered_distances = sort_and_filter_distances(distances)
          smoothed_distances = smooth_data(filtered_distances)
          slope, second_derivative = calculate_derivatives(smoothed_distances)
          curvature_peaks = find_curvature_peaks(second_derivative)
          max_curvature_index, max_curvature_point = find_max_curvature_point(filtered_distances, second_derivative)
          epsilon = max_curvature_point
          neighb = NearestNeighbors(n_neighbors=6)
          nbrs=neighb.fit(i)
          distances,indices=nbrs.kneighbors(i)
          epsilon_list.append(epsilon)
          
          dbscan = DBSCAN(eps = epsilon, min_samples = 6).fit(i) # fitting the model
          print(f'epsilon:{epsilon}')
          labels = dbscan.labels_
          x = i['x_pos'].astype(float)
          y = i['y_pos'].astype(float)
          z = i['z_pos'].astype(float)
          
          # Number of clusters in labels, ignoring noise if present.
          n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
          n_noise_ = list(labels).count(-1)
          df2 = pd.DataFrame()
          particle_ID = list(range(0, b-1))
          df2['particle_id'] = np.array(particle_ID)
          df2['labels'] = np.array(labels)
          df2['x_pos']= np.array(x)
          df2['y_pos']= np.array(y)
          df2['z_pos']= np.array(z)


          #separating particles which are outside the cluster
          df3 = df2[df2['labels'] == -1]
          bond_ = []          
          #saprate df for monomers outside the clusters
          df4 = df3.query("0 <= particle_id < 500")
          for i in df4['particle_id']:
              for j in df4['particle_id']:
                  if i == j + 1:
                     
                      x_i = df4.loc[df4['particle_id'] == i, 'x_pos'].iloc[0]
                      y_i = df4.loc[df4['particle_id'] == i, 'y_pos'].iloc[0]
                      z_i = df4.loc[df4['particle_id'] == i, 'z_pos'].iloc[0]
            
                      x_j = df4.loc[df4['particle_id'] == j, 'x_pos'].iloc[0]
                      y_j = df4.loc[df4['particle_id'] == j, 'y_pos'].iloc[0]
                      z_j = df4.loc[df4['particle_id'] == j, 'z_pos'].iloc[0]
            
                      diff_x = x_j - x_i
                      diff_y = y_j - y_i
                      diff_z = z_j - z_i
            
                      diff1_x = diff_x * diff_x
                      diff1_y = diff_y * diff_y
                      diff1_z = diff_z * diff_z    
                      magnitude = math.sqrt(diff1_x + diff1_y + diff1_z)    
                      bond_.append( magnitude)
                      
                      
                      mean__ = np.mean(bond_)
                     
                      
                      
          print(len(bond_))           
          mean_bond_length.append(mean__)            
                      
                      
l = np.mean(mean_bond_length)        

#average bond lengths calculated for the polymer in absence of proteins for different Re
#average bond length for Re = 100 is 1.09055
#average bond length for Re = 200 is 1.09059
#average bond length for Re = 250 is 1.09072
#Average bond length for Re = 300 is 1.09093
#Average bond length for Re = 350 is 1.09126
#average bond length for Re = 400 is 1.09133
#average bond length for Re = 500 is 1.09795

#constants 
l_0 = 1.09093 #mean bond length for Re=0.6
K = 100


delta_l = [item - l_0 for item in mean_bond_length]# for delta_length
tension = [item * K for item in delta_l] #K = 100

force = [] #converting in pN 

for i in tension:
             Fi = i*4.1/3.4
             force.append(Fi)

force = force[-3000:] #taking last 3000 values from each simulation

print(np.mean(force)) 
print(np.std(force))          

#saving force in pN
df8 = pd.DataFrame(force)
df8.to_csv('/home/rohitks/co-condensation/github/test/test_force.dat')

# Plot force vs time plot
plt.figure(figsize=(12, 6))
plt.plot(force, label="Force", color="Purple")
plt.title("Force vs time (Conc: 84.5, Re: 0.6, Replicate: 1)")
plt.xlabel("Time")
plt.ylabel("Force")
plt.grid(True)
plt.savefig('/home/rohitks/co-condensation/github/test/homo_0.6_rep_1_force_vs_time.svg', format='svg')
plt.close()



