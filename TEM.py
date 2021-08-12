from numpy import *
import matplotlib.pyplot as plt
import time


#################### Modify these values #######################################
save_file_ref_number = 5
''' We will make a matrix representing varying refractive indices. The indices will vary in both axes.
We will first define the size of the upper, mid, and lower slabs. Then add in the left and right
walls by modifying the outter values of the middle slab.
'''

#First we define the slab thicknesses in the verticle direction
d_top = 5E-6
d_center_vertical = 4E-6
d_bottom = 1E-6

#We then define the width of the middle slab. d_left, d_right are the left and right wall thicknesses, respectively.
d_center_horiz = 6E-6
d_left = 2E-6
d_right = 2E-6

#Resolution value, res, is a scaling factor of how many entries will be in the matrix of our waveguide
res = 20

#Light wavelength
lambda_0 = 13E-6
k_0 = (2*pi)/lambda_0 # Incident Wave Number

#refractive index for each layer
n_top = 1
n_mid = 1.5 
n_bottom = 1
n_left = 1
n_right = 1

# need to normalize thickness of dielectrics to get ratio for calculating number of points
d_min = min(d_top, d_center_vertical, d_bottom, d_left, d_right, d_center_horiz)

#number of points for each layer
layer_top = (d_top/d_min) * res
layer_mid = (d_center_vertical/d_min) * res
layer_bottom = (d_bottom/d_min) * res
layer_left = (d_left/d_min) * res
layer_center = (d_center_horiz/d_min) * res
layer_right = (d_right/d_min) * res

#Number of points in waveguid, not including the boundary. 
L = int(layer_top + layer_mid + layer_bottom)
W = int(layer_left + layer_center + layer_right)

#Create new variables to be consistent with the convention of M rows x N columns
M = L #Vertical gridpoints
N = W #Horizontal gridpoints

#Refractive index array, contains refractive index for each of the L x W points (L rows, W columns)
n = zeros([L, W])

#Step size
dx = d_min/res #Take the smallest thickness and divide by resolution

#Initialize array n
for k in range(0, L):
    if k < layer_top:
        for i in range (0, W):
            n[k, i] = n_top
    elif k < layer_top + layer_mid:
        for i in range (0, W):
            n[k,i] = n_mid
    else:
        for i in range (0, W):
            n[k, i] = n_bottom
print "Refractive index array, n, initialized vertically" 

for k in range(int(layer_top), int(layer_top+layer_mid) - 1):
    for i in range(0, W):
        if i < d_left:
            n[k, i] = n_left
        elif i < layer_left + layer_center - 1:    
            0
        else:
            n[k, i] = n_right

eigenMatrix_dimension = M*(N-2)
eigenMatrix = zeros([M*(N-2), M*(N-2)]) #the M(N-2) x M(N-2) matrix we need to solve

#first create matrix "alpha" from matrix n, refer to documentation 
alpha = zeros([L,W])
for i in range (0,L):
    for j in range(0,W):
        alpha[i,j] = (n[i,j]*dx*k_0)**2 - 4
        
#Initialize eigenMatrix diagonals
dummy_index_k = 0 #for indexing diagonals from (0,0) to (m(n-2), m(n-2))
for i in range (2 - 1,(N-1) - 1):
    for j in range(1 - 1,M - 1): #these range limits are determined by equations 4 and 5 in the documentation
        eigenMatrix[dummy_index_k, dummy_index_k] = alpha[i,j]
        dummy_index_k += 1
        
'''
#Use numpy to solve eigenvalues and eigenvectors
eig_vals = linalg.eig(A)[0]
eig_vects = linalg.eig(A)[1]
n_effective = sqrt(eig_vals)/(dx*k_0) #effective index

#Arrays for exporting to excel
eig_vects_export = eig_vects
index_bound = len(eig_vals)
k = 0
n_max = max(n_top, n_bottom)
while (k < index_bound):
    if eig_vals[k] < 0 or (n_effective[k] < n_max):
        print "Deleting %.3f, n_eff=%.3f" % (eig_vals[k], n_effective[k])
        eig_vects_export = delete(eig_vects_export, k, 1)
        eig_vals = delete(eig_vals, k)
        n_effective = delete(n_effective, k)
        index_bound = len(eig_vals)
    else:
        k += 1

#Code Below plots on different figures
x_axis = linspace(0, (L-1)*dx, L)

for k in range(0, len(eig_vals)):
    if eig_vals[k] < 0:
        continue
    plt.figure()
    plt.plot(x_axis, eig_vects_export[:,k])
    plt.title(r'$n_{eff}$ = %.3e $\lambda$ = %.3e' % (n_effective[k], lambda_0))
    plt.legend(loc='upper right')
    plt.pause(.001)
    plt.draw()


#The code below is a plot for on the same figurer
plt.figure()
for k in range(0, len(eig_vals)):
    if eig_vals[k] < 0:
        continue
    #plt.clf()
    plt.plot(x_axis, eig_vects_export[:,k], label="%.2e" % eig_vals[k])
    plt.title(r'$d_{mid}$ = %.3e $\lambda$ = %.3e' % (d_center_vertical, lambda_0))
    plt.legend(loc='upper right')
    plt.pause(1)
    plt.draw()


#Saving a csv of data

folder = 'run_%d' % save_file_ref_number
save_string_1 = "n_effective_%d" % save_file_ref_number
save_string_2 = "e_vects%d" % save_file_ref_number
savetxt(r"C:/Users/Panya/Documents/PHY409/Final Project/%s/%s.csv" % (folder, save_string_1), n_effective, fmt='%.3e', delimiter=',')
savetxt(r"C:/Users/Panya/Documents/PHY409/Final Project/%s/%s.csv" % (folder, save_string_2), eig_vects_export, fmt='%.3e', delimiter=',')
'''