from numpy import *
import matplotlib.pyplot as plt

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
plt.savefig('./Images/kvsPk.png')
plt.show()
