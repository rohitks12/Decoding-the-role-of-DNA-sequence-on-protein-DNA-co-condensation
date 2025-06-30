##################################################################################################################################################################
#In the following code, we are calculatin the probability of each monomer inside the condensate.
#To do this, we first identify particles belonging to the condensate using DBSCAN alogrithm.
#Next, we assign 0 to monomers outside the condensate, and 1 to monomers inside the condensate.
#Finally we calculate the probability of the monommer being inside the condensate by averaging over the scores computed from 3000 configuration files and store the data.
#This probability is used in the code for interfacial affinity analysis.
####################################################################################################################################################################

#importing libraries
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
from sklearn.neighbors import NearestNeighbors
from scipy.signal import savgol_filter
from natsort import natsorted
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


matrix = []
epsilon_list = []


# assigning directory
directory = '/home/rohitks/co-condensation/github/output_files/'


for filename in natsorted(os.listdir(directory)):
    f = os.path.join(directory, filename)

    if os.path.isfile(f):
        i = pd.read_csv(f, sep=" ", header=None, names=["index", "x_pos", "y_pos", "z_pos", "type"])
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
        nbrs = neighb.fit(i) 
        distances, indices = nbrs.kneighbors(i)
        epsilon_list.append(epsilon)

        dbscan = DBSCAN(eps=epsilon, min_samples=6).fit(i)  
        print(f'epsilon:{epsilon}')
        labels = dbscan.labels_
        x = i['x_pos'].astype(float)
        y = i['y_pos'].astype(float)
        z = i['z_pos'].astype(float)
        n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
        n_noise_ = list(labels).count(-1)
        df2 = pd.DataFrame()
        particle_ID = list(range(0, b - 1))
        df2['particle_id'] = np.array(particle_ID)
        df2['labels'] = np.array(labels)
        df2['x_pos'] = np.array(x)
        df2['y_pos'] = np.array(y)
        df2['z_pos'] = np.array(z)

        df3 = df2.query("-1 < particle_id < 500")
        labels__ = df3['labels']
        result = [1 if num > -1 else 0 for num in labels__]
        matrix.append(result)

##########################################################################################################################################################

df4 = pd.DataFrame(matrix)
df5 = df4.tail(3000)

probability = df5.mean(axis=0)
probability.to_csv('/home/rohitks/co-condensation/github/test/prob_0.002_rep_1.dat')

plt.figure(figsize=(12, 6))
plt.plot(probability, label="Probability", color="blue")
plt.title("Probability Array with Transitions (Conc: 0.002, Re: 0.6, Replicate: 1)")
plt.xlabel("Index")
plt.ylabel("Probability")
plt.grid(True)
plt.savefig('/home/rohitks/co-condensation/github/test/prob_profile.svg', format='svg')
plt.close()

##############################################################################################################################################################################################################

