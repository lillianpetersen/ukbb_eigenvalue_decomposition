## To install google cloud:
```
curl https://sdk.cloud.google.com | bash  
gcloud components install gsutil  
gcloud auth login
```

in ~/.bashrc:

```
EXPORT PATH=$PATH: /Users/.../google-cloud-sdk/bin
```

#### To list items in a bucket: 
```
gsutil ls (bucket)
```


## Running a vm:

### gcloud version
##### Make the Machine
``` 
gcloud compute instances create NAME --machine-type n1-standard-8  
gcloud compute ssh NAME
``` 

##### Install Hail on the Linux machine:
```
https://hail.is/docs/0.2/install/linux.html
```

##### Stop the Machine
```
gcloud compute instances stop NAME	(stops machine, but stores data)
```

##### Delete the Machine
```
gcloud compute instances delete NAME	(deletes entire machine)  
```

### dataproc version
``` 
hailctl dataproc start lkp --region us-central1 --requester-pays-allow-annotation-db --requester-pays-allow-all --pkgs matplotlib --num-workers 2 --num-preemptible-workers 4

hailctl dataproc modify lkp --region us-central1 --[how to modify]

hailctl dataproc connect --zone us-central1-a lkp notebook

hailctl dataproc list --region us-central1

hailctl dataproc stop lkp --region us-central1
```
### Megamem
$42/hour: careful!!

```
gcloud compute instances create lkpm --machine-type m1-ultramem-160
gcloud config set compute/zone us-central1-a
gcloud compute ssh lkpm
pip install -U hail
gcloud compute instances stop
gcloud compute instances delete
```


## In Python
#### Get correct java for Hail
```
export JAVA_HOME=/Library/Java/JavaVirtualMachines/jdk1.8.0_261.jdk/Contents/Home/
```

#### Export all autosomal chromosomes
```
CHROMOSOMES = [str(i) for i in range(1, 23)]  
CHROMOSOME_GLOB = '{' + ','.join(CHROMOSOMES) + '}'

mt = hl.import_bgen(  
	f'gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183/imputed/ukb_imp_chr{CHROMOSOME_GLOB}_v3.bgen',  
	entry_fields=['GT'],  
	sample_file='gs://ukb31063/ukb31063.autosomes.sample')
```


#### Hail Matrix table to numpy array:
```
np_gt = np.array(mt.GT.n_alt_alleles().collect()).reshape(100,100)

np.array(mt.GT.n_alt_alleles().collect()).reshape(mt.count_rows(),mt.count_cols())
```

## Neale Lab Cloud Documentation:  
```
https://github.com/danking/hail-cloud-docs
```