#! /home/pandunhuang/anaconda3/envs/mappy/bin/python
import os
import argparse
import mappy as mp
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('--r1',required=True,help="R1 fastq")
parser.add_argument('--r2',required=True,help="R2 fastq")
parser.add_argument('--name',required=True)
parser.add_argument('--output',default="./data")


args = parser.parse_args()

R1 = args.r1
R2 = args.r2
name = args.name
output = Path(args.output)
output.mkdir(exist_ok=True) 

outputr1 = output/f"{name}_S1_L001_R1_001.fastq"
outputr2 = output/f"{name}_S1_L001_R2_001.fastq.gz"

with open(outputr1,"w",encoding='utf-8') as out: 
    for ID,Seq,Q,Tag in mp.fastx_read(R1,read_comment=True):
        barcodex,Q1 = Seq[:8],Q[:8]
        barcodey,Q2 = Seq[38:46],Q[38:46]
        umi,Q3 = Seq[46:58],Q[46:58]
        out.write(f"{ID} {Tag}\n{barcodex}{barcodey}{umi}\n+\n{Q1}{Q2}{Q3}\n")

R2 = Path(R2).absolute()
os.system(f"pigz {outputr1}")
os.system(f"ln -s {R2} {outputr2}")
