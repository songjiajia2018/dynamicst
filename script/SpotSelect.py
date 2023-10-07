import json
import os
import pysam
import argparse
import scanpy as sc
import numpy as np
import pandas as pd
from pathlib import Path
from collections import Counter
from joblib import Parallel,delayed
from scipy.sparse import csr_matrix


parser = argparse.ArgumentParser()
parser.add_argument('--bam',required=True,help="bam file")
parser.add_argument('--barcode',default=None,help="spot list file under tissue")
parser.add_argument('--samplename',required=True,help="sample name")
parser.add_argument('--analysisdir',required=True,help="Analysis Dir")

args = parser.parse_args()

AnalysisDir             = Path(args.analysisdir)
SampleName              = args.samplename
Bamfile                 = args.bam
summary                 = AnalysisDir/SampleName/"Gene/Summary.csv"
barcodefile             = args.barcode
configfile              = AnalysisDir/"config.json"
featuresfile            = AnalysisDir/SampleName/"Gene/filtered/features.tsv"
bn = dict(Counter([i.strip() for i in open(barcodefile,encoding='utf-8')]))
barcodelist = [b.replace("_","") for b in bn if bn[b]==1]

Summary = pd.read_csv(summary,index_col=0,header=None).to_dict()[1]

def BarcodeInfo(Bamfile):
    d = {}
    barcodes_list = []
    with pysam.AlignmentFile(Bamfile, "rb") as samfile:
        for read in samfile:
            barcode = read.get_tag("CB").replace("_","")
            geneid  = read.get_tag("GX")
            umi     = read.get_tag("UB")
            if geneid != "-" and barcode != "-":
                if barcode in barcodelist:
                    barcodes_list.append([barcode,umi,geneid])
                d.setdefault(barcode,0)
                d[barcode] += 1
    return d,barcodes_list

def downsample(frac, df):
    """
    对DataFrame进行下采样，并计算UMI饱和度、spot的基因中位数

    参数：
    frac：下采样比例，介于0和1之间
    df：包含barcode、umi和geneid等相关信息的DataFrame

    返回值：
    umis_saturation：UMI饱和度
    avg_reads_per_spot：每个spot的平均读取数
    median_gene_per_spot：spot基因的中位数
    """

    # 如果下采样比例为0，则直接返回0
    if frac == 0:
        return 0, 0, 0
    d = {}
    n_deduped_reads = 0
    N_umis = 0
    N_reads = 0

    # 对DataFrame进行下采样
    dt = df.sample(frac=frac, random_state=0)
    median_gene = {}
    spots = set(dt[0].values)
    total_reads = frac * NumberOfReads

    # 遍历下采样后的数据
    for value in dt.values:
        bc, umi, gene = value
        median_gene.setdefault(bc, {})
        median_gene[bc][gene] = 1
        key = "".join(value)
        d.setdefault(key, 0)
        d[key] += 1

    # 计算UMI饱和度、读取饱和度和每个位置基因数量的中位数
    for key in d:
        if d[key] == 1:
            n_deduped_reads += 1
        N_umis += 1
        N_reads += d[key]

    umis_saturation = (1 - n_deduped_reads / N_umis)
    reads_saturation = (1 - n_deduped_reads / N_reads)
    median_gene2 = {bc: len(median_gene[bc]) for bc in median_gene}
    Median_Gene_Per_Spot = pd.DataFrame(median_gene2.items())[1].median()
    
    return umis_saturation, total_reads / len(spots), int(Median_Gene_Per_Spot)

def FilterMatrix():
    with open(str(AnalysisDir/SampleName/"Gene/filtered/barcodes.tsv"),"w",encoding='utf-8') as out:
        for barcode in barcodes:
            out.write(f"{barcode}\n")
    matrix = open(str(AnalysisDir/SampleName/"Gene/filtered/matrix.mtx"),"w")
    matrix.write("%%MatrixMarket matrix coordinate integer general\n")
    matrix.write("%\n")
    os.system(f"pigz {AnalysisDir/SampleName/'Gene/raw/*'}")
    adata = sc.read_10x_mtx(str(AnalysisDir/SampleName/"Gene/raw"))
    adata.obs.index = adata.obs.index.str.replace("_","")
    df_all_gene = adata.to_df()
    df_all_gene = df_all_gene.astype(int)
    df_all_gene = df_all_gene.loc[barcodes,features]
    sparse_matrix = csr_matrix(df_all_gene.values)
    nonzero_elements = sparse_matrix.nonzero()
    row_indices = nonzero_elements[0]+1  # 非零元素的行索引
    column_indices = nonzero_elements[1]+1  # 非零元素的列索引
    values = sparse_matrix.data  # 非零元素的值
    matrix.write(f"{len(features)} {len(barcodes)} {len(values)}\n")
    for gindex,bindex,value in zip(column_indices,row_indices,values):
        matrix.write(f"{gindex} {bindex} {value}\n")
    matrix.close()
    return df_all_gene
    

d,barcodes_list = BarcodeInfo(Bamfile)
df_chip         = pd.DataFrame(d.items())
df_chip.columns = ['barcode','count']
df_chip         = df_chip.set_index("barcode")
df_barcodes     = pd.DataFrame(barcodes_list)


NumberOfReads = int(Summary['Number of Reads'])
ValidBarcodes = Summary['Reads With Valid Barcodes']
downsample_saturation = Parallel(n_jobs=10)(delayed(downsample)(frac,df_barcodes) for frac in [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
df_saturation = pd.DataFrame(downsample_saturation,columns=['Saturation','Mean_Reads_Per_Spot','Median_Gene_Per_Spot'])

Saturation = df_saturation['Saturation'].values[-1]
metrics = {}

###
df_under_tissue = df_chip.loc[list(set(barcodelist)&set(df_chip.index))]

###Sequencing
metrics['Number of Reads']                      = NumberOfReads
metrics['Valid Barcodes']                       = ValidBarcodes
metrics['Sequencing Saturation']                = Saturation

barcodes = barcodelist
features = [i.split("\t")[1] for i in open(featuresfile)]
df_all_gene = FilterMatrix()
dt = df_all_gene.apply(lambda e:sum([1 for i in e if i !=0]),axis=1)
###Spots
metrics['Fraction Reads in Spots Under Tissue'] = f"{df_under_tissue['count'].sum()/df_chip['count'].sum():.2%}"
metrics['Number of Spots Under Tissue']         = len(barcodelist)
metrics['Mean Reads per Spot']                  = int(NumberOfReads/len(barcodelist))
metrics['UMIs in Spots']                        = int(df_all_gene.sum().sum())
metrics['Mean Reads Under Tissue per Spot']     = int(df_under_tissue['count'].mean())
metrics["Median UMI Counts per Spot"]           = int(df_all_gene.sum(axis=1).median())
metrics["Mean UMI Counts per Spot"]             = int(df_all_gene.sum(axis=1).mean())
metrics['Median Genes per Spot']                = int(dt.median())
metrics['Mean Genes per Spot']                  = int(dt.mean())
metrics['Total Genes Detected']                 = (df_all_gene.sum(axis=0) != 0 ).sum()

Summary.update(metrics)
pd.DataFrame(Summary.items()).to_csv(AnalysisDir/SampleName/"Gene/Summary.csv",index=False,header=False)
####修正中位基因数
diff_mgene = int(dt.median())-df_saturation['Median_Gene_Per_Spot'].values[-1]
df_saturation['Median_Gene_Per_Spot'] = df_saturation['Median_Gene_Per_Spot'].apply(lambda e:e+diff_mgene if e!=0 else e)
df_saturation.to_csv(AnalysisDir/f"{SampleName}_saturation.csv")
