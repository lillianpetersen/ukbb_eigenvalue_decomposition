#!/usr/bin/env python
# coding: utf-8

import hail as hl

from hail.plot import show
from os import path
import os
hl.plot.output_notebook()

import scipy.linalg
import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd
import time

# Set parameters
pct = 20
nsnps = 142909
nindv = 72174
print('\n'+str(pct)+'pct: '+str(nsnps)+' snps, '+str(nindv)+' indv\n')
    
imeta = 0

# Read grm
print('\nREADING GRM\n')
grm_file = "gs://ukb-gt/"+str(pct)+"pct/grm_"+str(nsnps)+"_x_"+str(nindv)+"_meta"+str(imeta)+".mt"
print(grm_file)
grm = hl.linalg.BlockMatrix.read(grm_file)
print('\nSUCCESSFULY READ GRM!\n')

# Convert to numpy
print('\nCONVERTING GRM TO NUMPY...\n')
np_grm = grm.to_numpy()
print('\nSUCCESSFULLY CONVERTED GRM TO NUMPY!\n')

# Compute eigenvalue decomposition
print('\nCOMPUTE EIGENVALUES...\n')
time_start = time.time()
eigenvals = scipy.linalg.eigvalsh(np_grm)
time_end = time.time()
eigenval_time = time_end - time_start
print(eigenval_time)
print('\nSUCCESSFULLY COMPUTED EIGENVALUES!\n')

# Save eigenvalues
np.save("eigenvals_"+str(pct)+"pct_"+str(nsnps)+"_x_"+str(nindv)+".npy", eigenvals)
np.save("/tmp/eigenvals_"+str(nsnps)+"_x_"+str(nindv)+".npy", eigenvals)
hl.hadoop_copy("file:///tmp/eigenvals_"+str(nsnps)+"_x_"+str(nindv)+".npy", "gs://ukb-gt/"+str(pct)+"pct/eigenvals_"+str(nsnps)+"_x_"+str(nindv)+".npy")
print('\nSAVED EIGENVALUES!\n')

# Plot
lmda = nindv / nsnps
lmdap = (1 + np.sqrt(lmda))**2
lmdam = (1 - np.sqrt(lmda))**2
x = np.arange(lmdam, lmdap, 0.001)
y = (1/(2*math.pi)) * np.sqrt((lmdap-x)*(x-lmdam)) / (lmda*x)

plt.clf()
plt.hist(eigenvals[1:], bins=1000, density = True)
plt.plot(x, y, '-b', label = 'Marchenko-Pastur Distrubution')
plt.title('Eigenvalues for '+str(nindv)+' Individuals and '+str(nsnps)+' SNPs')
plt.xlabel('Eigenvalues')
plt.ylabel('Density')
plt.xlim([0,5])
plt.grid(True)
plt.legend()
plt.savefig('eigenval_distribution_'+str(pct)+'pct_'+str(nsnps)+'_x_'+str(nindv)+'.pdf')
plt.savefig('/tmp/eigenval_dist_'+str(nindv)+'_x_'+str(nsnps)+'.pdf')
hl.hadoop_copy('file:///tmp/eigenval_dist_'+str(nindv)+'_x_'+str(nsnps)+'.pdf', 'gs://ukb-gt/'+str(pct)+'pct/eigenval_dist_'+str(nindv)+'_x_'+str(nsnps)+'.pdf')
print('\nCREATED PLOT!\n')

