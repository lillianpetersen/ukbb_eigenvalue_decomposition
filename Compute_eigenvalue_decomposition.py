#!/usr/bin/env python
# coding: utf-8

# In[84]:


import hail as hl

from hail.plot import show
from os import path
import os
import pandas as pd
import scipy.linalg
import numpy as np

hl.plot.output_notebook()


# In[85]:


#Useful command to specify nuimber of partitions:
def read_with_partitions(path, n_parts):
     mt = hl.read_matrix_table(path)
     return hl.read_matrix_table(path, _intervals=mt._calculate_new_partitions(n_parts))


# In[86]:


# Set parameters
pct = 20
print(str(pct)+"pct")


# In[87]:


#Specify file names and locations:
genomt_file = "gs://ukb31063/ukb31063.genotype.mt"

withdrawn_file = "gs://ukb31063/ukb31063.withdrawn_samples.ht"
gwassamples_file ="gs://ukb31063/ukb31063.neale_gwas_samples.both_sexes.ht"
sampleqc_file = "gs://ukb31063/ukb31063.sample_qc.tsv.bgz"
snpqc_file = "gs://ukb-gt/ukb_snp_qc_forHail_chr1to22.txt"

gwasvariants_file = "gs://ukb31063/ukb31063.neale_gwas_variants.ht"

phenos_file = "gs://ukb31063/ukb31063.PHESANT_January_2019.both_sexes.tsv.bgz"
assessment_file = "gs://ukb-gt/ukb_acentre.csv"
covariates_file = "gs://ukb31063/ukb31063.neale_gwas_covariates.both_sexes.tsv.bgz"

pcvar_file = "gs://nbaya/sex_linreg/ukb31063.*_phenotypes.both_sexes.reg1.tsv.bgz"


# In[88]:


# Define phenotype and assessment centre classification:
withdrawn = hl.read_table(withdrawn_file)
gwassamples = hl.read_table(gwassamples_file)


# In[89]:


# Import table of SNPs used in PCA conducted by the UKB:
snpqc = hl.import_table(snpqc_file, delimiter=" ", impute=True, key='rs_id')
snpqc = snpqc.filter(snpqc.in_PCA == 1, keep=True )  # only keep variants that have in_PCA == 1 (already passed LD pruning)
snpqc = snpqc.annotate(locus=hl.locus(hl.str(snpqc.chromosome), snpqc.position, reference_genome='GRCh37'))
snpqc = snpqc.key_by(snpqc.locus)


# In[90]:


# Read ukbb genotype MatrixTable and filter
ukbb = hl.read_matrix_table(genomt_file)
ukbb = ukbb.filter_cols(hl.is_defined(gwassamples[ukbb.col_key]), keep=True) #Keep only Neale lab GWAS samples
ukbb = ukbb.filter_cols(hl.is_defined(withdrawn[ukbb.col_key]), keep=False) #withdrawn samples were already removed actually
ukbb = ukbb.sample_cols(pct/100.) # (Lillian) filter number of individuals
ukbb = hl.variant_qc(ukbb)
ukbb = ukbb.filter_rows( (hl.min(ukbb.variant_qc.AF[0], ukbb.variant_qc.AF[1]) > 0.01) & ( ukbb.variant_qc.call_rate > 0.95), keep=True) # MAF > 0.01 and SNP call rate at least 95%
ukbb = hl.sample_qc(ukbb)
ukbb = ukbb.filter_cols( ukbb.sample_qc.call_rate > 0.98, keep=True) # sample missingness must be < 0.02
ukbb = hl.variant_qc(ukbb)
ukbb = ukbb.filter_rows( ukbb.variant_qc.call_rate > 0.98 , keep=True) # SNP missingness must be < 0.02 and hwe pvalue must be < 10^-6
ukbb = ukbb.filter_rows( ukbb.variant_qc.p_value_hwe < 0.0000000001 , keep=False) # SNP missingness must be < 0.02 and hwe pvalue must be < 10^-6
ukbb.count()


# In[91]:


# Save Metadata
metadata = pd.DataFrame({'nindv': [ukbb.count()[0]], 'nsnps': [ukbb.count()[1]], 'processing': ['SNPs from Neale Lab GWAS study']})

# Write UKB MT to file:
imeta = 0
while(True):
    if not path.exists("gs://ukb-gt/"+str(pct)+"pct/ukbQC_"+str(ukbb.count()[0])+"_x_"+str(ukbb.count()[1])+"_metadata"+str(imeta)+".txt"):
        metadata.to_csv("gs://ukb-gt/"+str(pct)+"pct/ukbQC_"+str(ukbb.count()[0])+"_x_"+str(ukbb.count()[1])+"_metadata"+str(imeta)+".txt", index=False)
        break
        
ukbQC_file = "gs://ukb-gt/"+str(pct)+"pct/ukbQC_"+str(ukbb.count()[0])+"_x_"+str(ukbb.count()[1])+"_meta"+str(imeta)+".mt"
ukbb.write(ukbQC_file, overwrite=True)
ukbb.count()


# In[92]:


#Read in UKB and with specified number of partitions:
ukbb = read_with_partitions(ukbQC_file, n_parts=500)
ukbb.count()


# In[93]:


# Filter SNPs to UKBB PCA SNPs
ukbb = ukbb.filter_rows(hl.is_defined(snpqc[ukbb.locus]))
ukbb.count()


# In[94]:


# Write the MTs of UKB and 1KG for PCA to file:
ukbPCA_file = "gs://ukb-gt/"+str(pct)+"pct/ukbPCA_"+str(ukbb.count()[0])+"_x_"+str(ukbb.count()[1])+"_meta"+str(imeta)+".mt"
ukbb.write(ukbPCA_file, overwrite=True)
print(ukbPCA_file)


# In[95]:


grm = hl.genetic_relatedness_matrix(ukbb.GT)
grm_file = "gs://ukb-gt/"+str(pct)+"pct/grm_"+str(ukbb.count()[0])+"_x_"+str(ukbb.count()[1])+"_meta"+str(imeta)+".mt"
grm.write(grm_file, overwrite=True)


# In[81]:


grm_np = grm.to_numpy()


# In[82]:


nsnps = ukbb.count()[0]
nindv = ukbb.count()[1]
eigenvals = scipy.linalg.eigvalsh(grm_np)
np.save("/tmp/eigenvals_"+str(nsnps)+"_x_"+str(nindv)+"_meta"+str(imeta)+".npy", eigenvals)
hl.hadoop_copy("file:///tmp/eigenvals_"+str(nsnps)+"_x_"+str(nindv)+"_meta"+str(imeta)+".npy", "gs://ukb-gt/"+str(pct)+"pct/eigenvals_"+str(nsnps)+"_x_"+str(nindv)+"_meta"+str(imeta)+".npy")


# In[83]:


import matplotlib.pyplot as plt
import math

print(nsnps, nindv)
lmda = nindv / nsnps
lmdap = (1 + np.sqrt(lmda))**2
lmdam = (1 - np.sqrt(lmda))**2
x = np.arange(lmdam, lmdap, 0.001)
y = (1/(2*math.pi)) * np.sqrt((lmdap-x)*(x-lmdam)) / (lmda*x)

plt.clf()
plt.hist(eigenvals[1:], bins=400, density = True)
plt.plot(x, y, '-b', label = 'Marchenko-Pastur Distrubution')
plt.title('Eigenvalues for '+str(nindv)+' Individuals and '+str(nsnps)+' SNPs')
plt.xlabel('Eigenvalues')
plt.ylabel('Density')
plt.xlim([0,2.25])
plt.grid(True)
plt.legend()
plt.savefig('/tmp/eigenval_dist_'+str(nindv)+'_x_'+str(nsnps)+'_meta'+str(imeta)+'.pdf')
hl.hadoop_copy('file:///tmp/eigenval_dist_'+str(nindv)+'_x_'+str(nsnps)+'_meta'+str(imeta)+'.pdf', 'gs://ukb-gt/'+str(pct)+'pct/eigenval_dist_'+str(nindv)+'_x_'+str(nsnps)+'_meta'+str(imeta)+'.pdf')
plt.show()


# In[ ]:




