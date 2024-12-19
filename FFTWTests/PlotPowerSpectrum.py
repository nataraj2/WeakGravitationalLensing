from numpy import *
import matplotlib.pyplot as plt
import os


data = loadtxt("k_vs_Pk.txt");
k = data[:,0]
Pk = data[:,1]
Pk_exact = data[:,2]


plt.figure(1)
plt.loglog(k, Pk, 'k', label='Computed')
plt.loglog(k, Pk_exact, 'ok', label='Specified')
plt.xlabel(r'$k$',fontsize=15)
plt.ylabel(r'$P(k)$',fontsize=15)
plt.legend()
# Define the directory name
directory = "Images"

# Check if the directory exists
if not os.path.exists(directory):
    # Create the directory
    os.makedirs(directory)

plt.savefig('./Images/kvsPk.png')
plt.show()
