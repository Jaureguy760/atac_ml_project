import pandas as pd
#Dict approach
import glob
import os
import sys
from pathlib import Path
#from typing import Any, Optional
import numpy as np
import pandas as pd
import torch
from torch.utils.data import Dataset, DataLoader
#import pytorch_lightning as pl
import h5py


#Load genome loader
sys.path.append('/iblm/netapp/home/jjaureguy/genome-loader/src')

sample_path = os.listdir('/iblm/netapp/data3/jjaureguy/PRJEB18997/10_genos/10_genos_fastq')

from encode_data import encode_sequence, encode_from_fasta, encode_from_h5
from get_encoded import get_encoded_haps
from load_data import load_vcf
from write_h5 import write_encoded_genome
from write_h5 import write_read_depth

test_fasta= '/iblm/netapp/data4/jjaureguy/ref_genome/hg38_ref/genome.fa'
out_dir = "/iblm/netapp/data4/jjaureguy/out_dir"
treatments_path = os.listdir('/iblm/netapp/data3/jjaureguy/PRJEB18997/10_genos/10_genos_fastq/HPSI0114i-eipl_1/out/treatments')

import glob
bam_paths = []
peak_paths = []
treatment_paths = []



column_names = [
    #'patient id', 
                'treatments','bam paths', 'peak paths']

d = {}

#filenames = [f for f in filenames if (f.startswith("ERR") and f.lower().endswith("1"))]
for i in (sample_path):
    #print([i][0]
    path = f"/iblm/netapp/data3/jjaureguy/PRJEB18997/10_genos/10_genos_fastq/{i}/out/treatments/*/merged_qc.bam"
    peak_path = f"/iblm/netapp/data3/jjaureguy/PRJEB18997/10_genos/10_genos_fastq/{i}/out/treatments/merged.bed"
    for filename in glob.glob(path):
        bam_paths.append(filename)
    for filename in glob.glob(peak_path):
        peak_paths.append(filename)
    for filename in (treatments_path):
        treatment_paths.append(filename)
    d[i] = pd.DataFrame(bam_paths)
    d[i].columns = ['bam paths']
    #d[i]['patient id'] = i
    d[i]['peak paths'] = pd.Series(peak_paths)
    d[i]['peak paths']= d[i]['peak paths'].fillna(method='ffill')
    d[i]['treatments']= pd.Series(treatment_paths)
    d[i]['treatments'] = d[i]['treatments'].str.replace('_mb', '')
    d[i] = d[i].reindex(columns=column_names)
    bam_paths.clear()
    peak_paths.clear()
    treatment_paths.clear()


    sample_path = os.listdir('/iblm/netapp/data3/jjaureguy/PRJEB18997/10_genos/19_genos_fastq')


    import glob
bam_paths = []
peak_paths = []
treatment_paths = []



column_names = [
    #'patient id', 
                'treatments','bam paths', 'peak paths']

d_add = {}

#filenames = [f for f in filenames if (f.startswith("ERR") and f.lower().endswith("1"))]
for i in (sample_path):
    #print([i][0]
    path = f"/iblm/netapp/data3/jjaureguy/PRJEB18997/10_genos/19_genos_fastq/{i}/out/treatments/*/merged_qc.bam"
    peak_path = f"/iblm/netapp/data3/jjaureguy/PRJEB18997/10_genos/19_genos_fastq/{i}/out/treatments/merged.bed"
    for filename in glob.glob(path):
        bam_paths.append(filename)
    for filename in glob.glob(peak_path):
        peak_paths.append(filename)
    for filename in (treatments_path):
        treatment_paths.append(filename)
    d_add[i] = pd.DataFrame(bam_paths)
    d_add[i].columns = ['bam paths']
    #d[i]['patient id'] = i
    d_add[i]['peak paths'] = pd.Series(peak_paths)
    d_add[i]['peak paths']= d_add[i]['peak paths'].fillna(method='ffill')
    d_add[i]['treatments']= pd.Series(treatment_paths)
    d_add[i]['treatments'] = d_add[i]['treatments'].str.replace('_mb', '')
    d_add[i] = d_add[i].reindex(columns=column_names)
    bam_paths.clear()
    peak_paths.clear()
    treatment_paths.clear()

df = pd.concat(d).reset_index().drop(columns='level_1').rename(columns={'level_0': 'patient_id'})
df_19 = pd.concat(d_add).reset_index().drop(columns='level_1').rename(columns={'level_0': 'patient_id'})

frames = [df, df_19]

result = pd.concat(frames)

from tqdm.notebook import tqdm

cur_offset = 0
bed_dfs = {}
for (patient_id, bam_path, peak_file), pdf in tqdm(result.groupby(['patient_id','bam paths', 'peak paths'])):
    
    bed_df = pd.read_csv(peak_file, delimiter='\t', 
                         header=None, names=['chrom', 'start', 'end'], 
                         dtype={'chrom':'category', 'start':int, 'end': int})
    
    bed_df['bam paths'] = bam_path
    bed_df['peak paths'] = peak_path
    bed_dfs[patient_id] = bed_df
dataset_df = pd.concat(bed_dfs, names=['patient_id']).reset_index()

del dataset_df['level_1']

#prev = dataset_df['bam paths'].unique()
x = result['bam paths'].unique()

#c = np.setxor1d(x, prev)

for fn in tqdm(x):
    fnp = Path(fn)
    outdir = fnp.parent
    write_read_depth(str(fnp), str(outdir))


