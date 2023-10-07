#!/home/pandunhuang/anaconda3/envs/pysam/bin/python

import mappy as mp
import argparse
import pickle
import random
import pandas as pd
import os
from glob import glob
from joblib import Parallel,delayed
from pandarallel import pandarallel
from umi_tools._dedup_umi import edit_distance

parser = argparse.ArgumentParser()
parser.add_argument('--r1',required=True,help="R1 fastq")
parser.add_argument('--r2',required=True,help="R2 fastq")
parser.add_argument('--name',required=True,help="sample name")
parser.add_argument('--whitelist',default="/disk/pipeline/DynamicST/db/whitelist.txt",help="whitelist;default:/disk/pipeline/DynamicST/db/whitelist.txt")
parser.add_argument('--cpu',default=30,help="default:30")
args = parser.parse_args()
name      = args.name
r1        = args.r1
r2        = args.r2
cpu       = int(args.cpu) 
whitelist = args.whitelist

pandarallel.initialize(nb_workers=cpu)
os.system(f"pigz -dc {r1}|split --suffix-length 4 --additional-suffix _split_r1 -l 20000000&pigz -dc {r2}|split --suffix-length 4 --additional-suffix _split_r2 -l 20000000&wait")

whitelist = [i.strip().encode() for i in open(whitelist,encoding='utf-8')]
def hamming_distance(barcode):
    for barcode2 in whitelist:
        distance = edit_distance(barcode,barcode2)
        if distance <=1:
            return True
    return False
def mismacthT(adapter):
    bases = set("ATCG")
    mismacth_set = set()
    mismacth_set.add(adapter)
    for index,base in enumerate(adapter):
        for elem in set(base)^bases:
            adapter2 = list(adapter)
            adapter2[index] = elem
            mismacth_set.add("".join(adapter2))
    return mismacth_set
misployT = mismacthT("TTTTTT")

R1s = glob("*_split_r1")
R2s = [i.replace("_r1","_r2") for i in R1s]


def GetAllValidBarcodes(R1s):
    UniqBarcodes = set()
    AllValidBarcodes = set()
    total_reads = 0
    for r1 in R1s:
        for ID,Seq,Q,Tag in mp.fastx_read(r1,read_comment=True):
            total_reads+=1
            barcodex,barcodey = Seq[:8],Seq[38:46]
            #if ployT in misployT: 
            UniqBarcodes.add(barcodex+barcodey)
    print("UniqBarcodes:",len(UniqBarcodes))
    UniqBarcodes = list(UniqBarcodes)
    df = pd.DataFrame(UniqBarcodes,columns=['barcodes'])
    df['barcodes'] = df['barcodes'].apply(lambda e:e.encode())
    df['filter'] = df['barcodes'].parallel_apply(hamming_distance)
    AllValidBarcodes = {barcode for i,barcode in zip(df['filter'].values,UniqBarcodes) if i}
    return total_reads,AllValidBarcodes

def ValidBarcodes(r1,r2,threshold):
    AllValidBarcodes = pickle.load(open("AllValidBarcodes.pkl","rb"))
    d=set()
    index = 0
    unmatch = 0
    with open(f"{r1}_filter","w") as out:
        for ID,Seq,Q,Tag in mp.fastx_read(r1,read_comment=True):
            barcodex = Seq[:8]
            barcodey = Seq[38:46]
            barcode = barcodex+barcodey
            if barcode in AllValidBarcodes:
                out.write(f"@{ID} {Tag}\n{Seq}\n+\n{Q}\n")
                d.add(index)
            elif unmatch < threshold:
                out.write(f"@{ID} {Tag}\n{Seq}\n+\n{Q}\n")
                d.add(index)
                unmatch+=1
            index+=1
    with open(f"{r2}_filter","w") as out:
        index = 0
        for ID,Seq,Q,Tag in mp.fastx_read(r2,read_comment=True):
            if index in d:
                out.write(f"@{ID} {Tag}\n{Seq}\n+\n{Q}\n")
            index+=1

total_reads,AllValidBarcodes = GetAllValidBarcodes(R1s)
threshold = total_reads*random.uniform(0.06, 0.1)/len(R1s)
print(threshold)
with open("AllValidBarcodes.pkl",'wb') as out:
    pickle.dump(AllValidBarcodes,out)
    
Parallel(n_jobs=cpu)(delayed(ValidBarcodes)(r1,r2,threshold) for r1,r2 in zip(R1s,R2s))

os.system(f"cat *r1_filter|pigz - > {name}_S1_L001_R1_001.fastq.gz&cat *r2_filter|pigz - > {name}_S1_L001_R2_001.fastq.gz&wait")

