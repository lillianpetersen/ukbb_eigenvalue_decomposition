import numpy as np
import matplotlib.pyplot as plt
import math

pct = 1
nsnps = 145071
nindv = 3679
#pct = 10
#nsnps = 144205
#nindv = 35946
#pct = 5
#nsnps = 145061
#nindv = 17799
imeta = 0

eigenvals = np.load(str(pct)+'pct/eigenvals_'+str(nsnps)+'_x_'+str(nindv)+'_meta'+str(imeta)+'.npy')

print(nsnps, nindv)
lmda = nindv / nsnps
lmdap = (1 + np.sqrt(lmda))**2
lmdam = (1 - np.sqrt(lmda))**2
x = np.arange(lmdam, lmdap, 0.001)
y = (1/(2*math.pi)) * np.sqrt((lmdap-x)*(x-lmdam)) / (lmda*x)

plt.clf()
plt.hist(eigenvals[1:], bins=int(nindv/35), density = True)
plt.plot(x, y, '-b', label = 'Marchenko-Pastur Distrubution')
plt.title('Eigenvalues for '+str(nindv)+' Individuals and '+str(nsnps)+' SNPs')
plt.xlabel('Eigenvalues')
plt.ylabel('Density')
plt.xlim([0,2.75])
plt.grid(True)
plt.legend()
plt.savefig(str(pct)+'pct/eigenval_dist_'+str(pct)+'pct_'+str(nsnps)+'_x_'+str(nindv)+'_meta'+str(imeta)+'.pdf')
