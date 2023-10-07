import numpy as np
import os
import seaborn as sns
import argparse
import cv2
import json
import h5py
import pandas as pd
from pathlib import Path
import scanpy as sc
import matplotlib.pyplot as plt
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument('--softhome',required=True)
parser.add_argument('--image',required=True,help="Corrected HE image")
parser.add_argument('--inputdir',required=True,help="Gene mtx dir")
parser.add_argument('--sample',required=True,help="Sample name")
parser.add_argument('--org',default="mm10",help='default mm10')
parser.add_argument('--STARLog',required=True)
parser.add_argument('--RnaSeqMetrics',required=True)
parser.add_argument('--sampleinfo',required=True)
parser.add_argument("--report_temp",required=True)
parser.add_argument("--output",required=True)
parser.add_argument("--fastp",required=True)
parser.add_argument("--MetricsSummary",required=True)

args = parser.parse_args()

def cell_and_sequencing_summary(STARSummary):
    summary = {"cell":{},"sequencing":{}}
    df = pd.read_csv(STARSummary,index_col=0,header=None)
    df = df.dropna()
    df.columns=['index']
    d = df.to_dict()['index']
    summary['sequencing']["Reads Mapped to Genome"] = f"{d['Reads Mapped to Genome: Unique+Multiple']:.2%}"
    summary['sequencing']['Number of Reads'] = f"{int(d['Number of Reads']):,}"
    summary['sequencing']['Valid Barcodes'] = f"{d['Reads With Valid Barcodes']:.2%}"
    summary['sequencing']['Sequencing Saturation'] = f"{d['Sequencing Saturation']:.2%}"
    if d.get('Q30 Bases in CB+UMI'):
        summary['sequencing']['CB+UMI'] = f"""
    <tr>
      <th>Q30 Bases in CB+UMI</th>
      <th>{d['Q30 Bases in CB+UMI']:.2%}</th>
    </tr>
    """
    else:
        summary['sequencing']['CB+UMI'] = ""
    summary['sequencing']['Q30 Bases in RNA Read'] = f"{d['Q30 Bases in RNA read']:.2%}"
    summary['cell']['Estimated Number of Cells'] = f"{int(d['Estimated Number of Cells']):,}"
    summary['cell']['Fraction Reads in Cells'] = f"{d['Fraction of Unique Reads in Cells']:.2%}"
    summary['cell']['Mean Reads per Cell'] = f"{int(int(d['Number of Reads'])/int(d['Estimated Number of Cells'])):,}"
    summary['cell']['Median UMI Counts per Cell'] = f"{int(d['Median UMI per Cell']):,}"
    summary['cell']['Median Genes per Cell'] = f"{int(d['Median Gene per Cell']):,}"
    summary['cell']['Total Genes Detected'] = f"{int(d['Total Gene Detected']):,}"
    return summary

def mapping_summary(STARLog, RnaSeqMetrics):
    d={}
    summary = {}
    with open(STARLog, 'r',encoding='utf-8') as fh:
        for line in fh:
            if 'Number of input reads' in line:
                summary['Number of input reads'] = int(
                    line.strip().split('\t')[-1])
            if 'Uniquely mapped reads number' in line:
                summary['Uniquely mapped reads number'] = int(
                    line.strip().split('\t')[-1])
            if 'Number of reads mapped to multiple loci' in line:
                summary['Number of reads mapped to multiple loci'] = int(
                    line.strip().split('\t')[-1])
    with open(RnaSeqMetrics, 'r',encoding='utf-8') as fh:
        while True:
            line = fh.readline().strip()
            if line.startswith('total alignments'):
                summary['total alignments'] = int(line.split()[-1].replace(
                    ',', ''))
            if line.startswith('reads aligned'):
                summary['reads aligned'] = int(line.split()[-1].replace(
                    ',', ''))
            if line.startswith('aligned to genes'):
                summary['aligned to genes'] = int(line.split()[-1].replace(
                    ',', ''))
            if line.startswith('no feature assigned'):
                summary['no feature assigned'] = int(line.split()[-1].replace(
                    ',', ''))
            if line.startswith('exonic'):
                summary['exonic'] = int(line.split()[-2].replace(',', ''))
            if line.startswith('intronic'):
                summary['intronic'] = int(line.split()[-2].replace(',', ''))
            if line.startswith('intergenic'):
                summary['intergenic'] = int(line.split()[-2].replace(',', ''))
                break
    #d['Reads Mapped to Genome'] = f"{summary['reads aligned']/summary['Number of input reads']:.2%}"
    intergenic = summary['intergenic']/summary['Number of input reads']
    intronic = summary['intronic']/summary['Number of input reads']
    exonic = summary['exonic']/summary['Number of input reads']
    d['Reads Mapped Confidently to Genome'] = f"{intergenic+intronic+exonic:.2%}"
    d['Reads Mapped Confidently to Intergenic Regions'] = f"{intergenic:.2%}"
    d['Reads Mapped Confidently to Intronic Regions'] = f"{intronic:.2%}"
    d['Reads Mapped Confidently to Exonic Regions'] = f"{exonic:.2%}"
    return d

def Sample_summary(config):
    summary = {}
    d = json.load(open(config))
    summary['Sample ID'] = d['sample']
    summary['Sample Description'] = d.get("Sample Description","")
    summary['Chemistry'] = "Single Cell 3' v1"
    summary['Include introns'] = "False"
    summary['Transcriptome'] = Path(d['transcriptome']).parent.name
    summary['Pipeline Version'] = "DynamicST-1.0.0"
    return summary

def fastp_qc(fastp):
    summary={}
    d = json.load(open(fastp))
    summary["gc_content"] = f"{d['summary']['before_filtering']['gc_content']:.2%}"
    summary["too_short_reads"] = f"{d['filtering_result']['too_short_reads']:,}"
    summary["too_many_N_reads"] = f"{d['filtering_result']['too_many_N_reads']:,}"
    summary['q20_rate'] = f"{d['summary']['before_filtering']['q20_rate']:.2%}"
    return summary

def Get_adata():
    adata = sc.read_10x_mtx(str(Path(inputdir)/"filtered"))
    sc.pp.filter_genes(adata, min_cells=1)
    sc.pp.filter_cells(adata, min_genes=1)
    return adata

def get_tissue(adata,coord):
    d= adata.obs.to_dict()['n_genes']
    tissue=[]
    for barcode,(x,y) in coord.iterrows():
        umi_count = d.get(barcode,0)
        if umi_count > 0:
            tissue.append(barcode)
    return tissue

def get_ticks(HE):
    img = cv2.imread(HE)
    stepx = img.shape[0]/51
    stepy = img.shape[1]/51
    xcoord = []
    ycoord = []
    for x in range(1,51):
        for y in range(1,51):
            xcoord.append(x*stepx)
            ycoord.append(y*stepy)
    return xcoord,ycoord,img

def Spots_plot(img,xcoord,ycoord,coord):
    (Path(output)/"src").mkdir(exist_ok=True)
    plt.rcParams['figure.dpi']=100
    h,w = img.shape[:2]
    scale = h/w
    plt.figure(figsize=(10,10*scale))
    plt.imshow(img)
    plt.scatter(xcoord,ycoord,c="white",s=100,marker='s',edgecolors="black",alpha=0.2)
    x1,y1 = [],[]
    for i in coord.values:
        y,x,in_tissue,pxl_col_in_fullres,pxl_row_in_fullres = i
        if in_tissue == 1:
            x1.append(pxl_row_in_fullres)
            y1.append(pxl_col_in_fullres)
    plt.scatter(x1,y1,c="red",marker='s',s=100,edgecolors="black",alpha=0.2)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(f"{output}/src/Spots.png")

    
def Generate_spatial_directory():
    coord = pd.read_csv(Path(softhome)/'template/barcode_coordinate.txt',sep='\t',index_col=0)
    adata = Get_adata()
    tissue = get_tissue(adata,coord)
    xcoord,ycoord,img = get_ticks(HE)
    coord['in_tissue'] = [0 if b not in tissue else 1 for b in coord.index]
    coord["pxl_col_in_fullres"] = ycoord
    coord["pxl_row_in_fullres"] = xcoord
    Path(f"{inputdir}/../outs/spatial").mkdir(exist_ok=True,parents=True)
    coord = coord[['in_tissue','y','x','pxl_col_in_fullres','pxl_row_in_fullres']]
    coord.to_csv(f"{inputdir}/../outs/spatial/tissue_positions_list.csv",header=None)
    Spots_plot(img,xcoord,ycoord,coord)
    spot_diameter = coord['pxl_col_in_fullres'].values[1] - coord['pxl_col_in_fullres'].values[0]
    d={"tissue_hires_scalef":2000/6000,
    "tissue_lowres_scalef":600/6000,
    "fiducial_diameter_fullres":spot_diameter+30,
    "spot_diameter_fullres":spot_diameter}
    with open(f"{inputdir}/../outs/spatial/scalefactors_json.json","w") as out:
        json.dump(d,out)
    hires_img = cv2.resize(img,(2000,2000))
    lowres_img = cv2.resize(img,(600,600))
    cv2.imwrite(f"{inputdir}/../outs/spatial/tissue_hires_image.png",hires_img)
    cv2.imwrite(f"{inputdir}/../outs/spatial/tissue_lowres_image.png",lowres_img)

def Generate_H5():
    f = h5py.File(f"{inputdir}/../outs/filtered_feature_bc_matrix.h5", "w")
    f.attrs['library_ids'] = np.array([f'{sample}'],dtype="|S27")
    matrix = f.create_group("matrix")
    barcode = pd.read_csv(f"{inputdir}/filtered/barcodes.tsv.gz",sep="\t",header=None)[0].values
    matrix.create_dataset("barcodes", data=barcode)
    df_mtx = pd.read_csv(f"{inputdir}/filtered/matrix.mtx.gz",sep="\s+",header=None,skiprows=3)
    df_mtx[0] = df_mtx[0] - 1

    data = df_mtx[2].values
    matrix.create_dataset("data", data=data)


    features = matrix.create_group("features")
    df_features = pd.read_csv(f"{inputdir}/filtered/features.tsv.gz",header=None,sep='\t')

    features.create_dataset("_all_tag_keys", data=[b'genome'])
    features.create_dataset("feature_type", data=df_features[2].values)
    features.create_dataset("genome", data=[genome]*len(df_features))
    features.create_dataset("id", data=df_features[0].values)
    features.create_dataset("name",data=df_features[1].values)
    indices = df_mtx[0].values
    matrix.create_dataset("indices", data=indices)
    indptr = np.array([0]+list(df_mtx[1].value_counts().sort_index().cumsum().values))
    matrix.create_dataset("indptr", data=indptr)
    shape = np.array([len(df_features),len(df_mtx[1].unique())])
    matrix.create_dataset("shape", data=shape)
    f.close()

RnaSeqMetrics = args.RnaSeqMetrics
sampleinfo = args.sampleinfo
STARLog = args.STARLog
Temp = args.report_temp
genome = args.org
###
output = Path(args.output)

os.system(f"cp -r {Path(Temp).parent}/'static' {output}")
fastp = args.fastp
softhome = args.softhome
HE = args.image
inputdir = Path(args.inputdir)
STARSummary = inputdir/"Summary.csv"
sample = args.sample
MetricsSummary = args.MetricsSummary
Path(f"{output}/spatial").mkdir(exist_ok=True,parents=True)

total_summary = {}
total_summary["qc"] = fastp_qc(fastp)
total_summary.update(cell_and_sequencing_summary(STARSummary))
total_summary['mapping'] = mapping_summary(STARLog, RnaSeqMetrics)
total_summary['qc']['GC Content'] = total_summary['qc']['gc_content']
total_summary['qc']['Q20 Bases in RNA Read'] = total_summary['qc']['q20_rate']
del total_summary['qc']['gc_content']
del total_summary['qc']['q20_rate']
total_summary['sample'] = Sample_summary(sampleinfo)


metrics_summary = []
ms = {}
for key in total_summary:
    ms.update(total_summary[key])
web_summary={}
for key in ms:
    value = ms[key]
    if "Cells" in key:
        key = key.replace("Cells",'Spots')
    elif "Cell" in key:
        key = key.replace("Cell",'Spot')
    web_summary[key] = value

for key in ["Number of Reads","Valid Barcodes","Sequencing Saturation","GC Content","Q20 Bases in RNA Read","Q30 Bases in RNA Read","Estimated Number of Spots","Fraction Reads in Spots","Mean Reads per Spot","Median UMI Counts per Spot","Median Genes per Spot","Total Genes Detected","Reads Mapped to Genome","Reads Mapped Confidently to Genome","Reads Mapped Confidently to Intergenic Regions","Reads Mapped Confidently to Intronic Regions","Reads Mapped Confidently to Exonic Regions"]:
    metrics_summary.append([key,str(web_summary[key]).replace(',','')])

pd.DataFrame(metrics_summary).T.to_csv(MetricsSummary,index=False,header=False)

f = open(Temp,encoding='utf-8').read()
for key in web_summary:
    f = f.replace(f"$${key}$$",str(web_summary[key]))
with open(str(Path(output)/"web_summary.html"),"w",encoding='utf-8') as out:
    out.write(f)

Generate_H5()
Generate_spatial_directory()