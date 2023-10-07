import cv2
import os
import h5py
import json
import argparse
import tifffile
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path
import matplotlib.pyplot as plt
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument('--AnalysisDir',help='Analysis Directory',required=True)
parser.add_argument('--Sample',help="Sample Name",required=True)
parser.add_argument('--SoftHome',help="Software Analysis Directory;default:/disk/pipeline/DynamicST",default="/disk/pipeline/DynamicST")
args = parser.parse_args()


def Spot_and_sequencing_summary(STARSummary):
    global NumberOfReads
    summary = {"Spot":{},"sequencing":{}}
    df = pd.read_csv(STARSummary,index_col=0,header=None)
    #df.index = df.index.str.replace("Cells","Spots").str.replace("Cell","Spot")
    df = df.dropna()
    df.columns=['index']
    d = df.to_dict()['index']
    NumberOfReads = int(d['Number of Reads'])
    summary['sequencing']["Reads Mapped to Genome"] = f"{float(d['Reads Mapped to Genome: Unique+Multiple']):.2%}"
    summary['sequencing']['Number of Reads']        = f"{NumberOfReads:,}"
    summary['sequencing']['Valid Barcodes']         = f"{float(d['Reads With Valid Barcodes']):.2%}"
    summary['sequencing']['Sequencing Saturation']  = f"{float(d['Sequencing Saturation']):.2%}"
    if d.get('Q30 Bases in CB+UMI'):
        summary['sequencing']['CB+UMI'] = f"""
    <tr>
      <th>Q30 Bases in CB+UMI</th>
      <th>{float(d['Q30 Bases in CB+UMI']):.2%}</th>
    </tr>
    """
    else:
        summary['sequencing']['CB+UMI'] = ""
    summary['sequencing']['Q30 Bases in RNA Read']          = f"{float(d['Q30 Bases in RNA read']):.2%}"
    summary['Spot']['Number of Spots Under Tissue']         = f"{int(d['Number of Spots Under Tissue']):,}"
    summary['Spot']['Fraction Reads in Spots Under Tissue'] = d['Fraction Reads in Spots Under Tissue']
    summary['Spot']['Mean Reads per Spot']                  = f"{int(float(d['Mean Reads per Spot'])):,}"
    summary['Spot']['Mean Reads Under Tissue per Spot']     = f"{int(float(d['Mean Reads Under Tissue per Spot'])):,}"
    summary['Spot']["UMIs in Spots"]                        = f"{int(float(d['UMIs in Spots'])):,}"
    summary['Spot']['Median UMI Counts per Spot']           = f"{int(float(d['Median UMI Counts per Spot'])):,}"
    summary['Spot']['Median Genes per Spot']                = f"{int(float(d['Median Genes per Spot'])):,}"
    summary['Spot']['Total Genes Detected']                 = f"{int(float(d['Total Genes Detected'])):,}"
    return summary

def mapping_summary(RnaSeqMetrics):
    d={}
    summary = {}
    with open(RnaSeqMetrics, 'r',encoding='utf-8') as fh:
        while True:
            line = fh.readline().strip()
            if line.startswith('exonic'):
                summary['exonic'] = int(line.split()[-2].replace(',', ''))
            if line.startswith('intronic'):
                summary['intronic'] = int(line.split()[-2].replace(',', ''))
            if line.startswith('intergenic'):
                summary['intergenic'] = int(line.split()[-2].replace(',', ''))
                break
    intergenic = summary['intergenic']/NumberOfReads
    intronic = summary['intronic']/NumberOfReads
    exonic = summary['exonic']/NumberOfReads
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
    summary['Chemistry'] = "Single Spot 3' v1"
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

def GetCoord(img,TotalSpot):
    img = cv2.imread(img)
    img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
    stepx = img.shape[0]/(TotalSpot+1)
    stepy = img.shape[1]/(TotalSpot+1)
    xcoord = []
    ycoord = []
    for x in range(1,(TotalSpot+1)):
        for y in range(1,(TotalSpot+1)):
            xcoord.append(x*stepx)
            ycoord.append(y*stepy)
    return img,xcoord,ycoord

def SpotsPlot(img,df):
    plt.figure(figsize=(14,14))
    plt.imshow(img)
    dt = df[df['in_tissue'] == 1]
    for y,x in dt[['pxl_col_in_fullres','pxl_row_in_fullres']].values:
        plt.scatter(x,y,color='red',s=40,marker='s',lw=0.4,edgecolors='black')
    plt.axis("off")
    Path(f"{AnalysisDir}/{sample}/outs/src").mkdir(exist_ok=True,parents=True)
    plt.savefig(f"{AnalysisDir}/{sample}/outs/src/Spots.png")

def MakeSpatial(img,barcodex,barcodey):
    try:
        barcode = pd.read_csv(str(AnalysisDir/f"{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"),sep='\t',header=None)
    except:
        barcode = pd.read_csv(str(AnalysisDir/f"{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv"),sep='\t',header=None)
    tissue = list(barcode[0].values)
    df = pd.DataFrame([[f"{bx.strip()}{by.strip()}",xindex,yindex] for xindex,bx in enumerate(open(barcodex,encoding='utf-8'),1) for yindex,by in enumerate(open(barcodey,encoding='utf-8'),1)])
    df.columns = ['barcode','x','y']
    df = df.set_index("barcode")
    #df.to_csv(f"{AnalysisDir}/barcode_coordinate.txt",sep='\t',encoding='utf-8')
    TotalSpot = df['x'].max()
    img,xcoord,ycoord = GetCoord(img,TotalSpot)
    #if TotalSpot == 96:
    #    df = pd.read_csv(str(SoftHome/'db/96barcode_coordinate.txt'),sep='\t',index_col=0)
    #else:
    #    df = pd.read_csv(str(SoftHome/'db/barcode_coordinate.txt'),sep='\t',index_col=0)

    df["pxl_col_in_fullres"] = ycoord
    df["pxl_row_in_fullres"] = xcoord
    df['in_tissue'] = [0 if b not in tissue else 1 for b in df.index]
    SpotsPlot(img,df)
    Path(f"{AnalysisDir}/{sample}/outs/spatial").mkdir(exist_ok=True,parents=True)
    df = df[['in_tissue','y','x','pxl_col_in_fullres','pxl_row_in_fullres']]
    df.to_csv(f"{AnalysisDir}/{sample}/outs/spatial/tissue_positions_list.csv",header=None)
    spot_diameter = df['pxl_col_in_fullres'].values[1] - df['pxl_col_in_fullres'].values[0]
    d={"tissue_hires_scalef":2000/6000,
    "tissue_lowres_scalef":600/6000,
    "fiducial_diameter_fullres":spot_diameter+30,
    "spot_diameter_fullres":spot_diameter}
    with open(f"{AnalysisDir}/{sample}/outs/spatial/scalefactors_json.json","w") as out:
        json.dump(d,out)
    hires_img = cv2.resize(img,(2000,2000))
    lowres_img = cv2.resize(img,(600,600))
    cv2.imwrite(f"{AnalysisDir}/{sample}/outs/spatial/tissue_hires_image.png",hires_img)
    cv2.imwrite(f"{AnalysisDir}/{sample}/outs/spatial/tissue_lowres_image.png",lowres_img)


def MakeH5File(genome):
    f = h5py.File(f"{AnalysisDir}/{sample}/outs/filtered_feature_bc_matrix.h5", "w")
    f.attrs['library_ids'] = np.array([f'{sample}'],dtype="|S27")
    matrix = f.create_group("matrix")
    try:
        barcode = pd.read_csv(f"{AnalysisDir}/{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",sep="\t",header=None)[0].values
    except:
        barcode = pd.read_csv(f"{AnalysisDir}/{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv",sep="\t",header=None)[0].values
    matrix.create_dataset("barcodes", data=barcode)
    try:
        df_mtx = pd.read_csv(f"{AnalysisDir}/{sample}/outs/filtered_feature_bc_matrix/matrix.mtx.gz",sep="\s+",header=None,skiprows=3)
    except:
        df_mtx = pd.read_csv(f"{AnalysisDir}/{sample}/outs/filtered_feature_bc_matrix/matrix.mtx",sep="\s+",header=None,skiprows=3)
    df_mtx[0] = df_mtx[0] - 1
    data = df_mtx[2].values
    matrix.create_dataset("data", data=data)
    features = matrix.create_group("features")
    try:
        df_features = pd.read_csv(f"{AnalysisDir}/{sample}/outs/filtered_feature_bc_matrix/features.tsv.gz",header=None,sep='\t')
    except:
        df_features = pd.read_csv(f"{AnalysisDir}/{sample}/outs/filtered_feature_bc_matrix/features.tsv",header=None,sep='\t')

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


AnalysisDir   = Path(args.AnalysisDir)
SoftHome      = Path(args.SoftHome)
sample        = args.Sample
img           = str(AnalysisDir/'image/HEadj.tif')
RnaSeqMetrics = str(AnalysisDir/"qual_summary/rnaseq_qc_results.txt")
saturation    = str(AnalysisDir/"saturation.csv")
sampleinfo    = str(AnalysisDir/"config.json")
STARSummary   = str(AnalysisDir/f"{sample}/Gene/Summary.csv")
output        = str(AnalysisDir/f"{sample}/outs/web_summary.html")
Temp          = str(SoftHome/'template/WPSQCReport_temp_v2.html')
fastp         = str(AnalysisDir/"fastp.json")
MetricsSummary = str(AnalysisDir/sample/'outs/metrics_summary.csv')

genome = Path(json.load(open(sampleinfo,encoding='utf-8'))['transcriptome']).parent.name
whitelist = json.load(open(sampleinfo,encoding='utf-8'))['whitelist']
if type(whitelist) == str:
    barcodex,barcodey = whitelist.split()
else:
    barcodex,barcodey = whitelist
MakeSpatial(img,barcodex,barcodey)
MakeH5File(genome)

total_summary = {}
total_summary["qc"] = fastp_qc(fastp)
total_summary.update(Spot_and_sequencing_summary(STARSummary))
total_summary['mapping'] = mapping_summary(RnaSeqMetrics)
total_summary['sample'] = Sample_summary(sampleinfo)

metrics_summary = []

ms = {}
for key in total_summary:
    ms.update(total_summary[key])
trans = {"GC Content":"gc_content","Q20 Bases in RNA Read":"q20_rate"}
print(ms)
for key in ["Number of Reads","Valid Barcodes","Sequencing Saturation","GC Content","Q20 Bases in RNA Read","Q30 Bases in RNA Read","Number of Spots Under Tissue","Fraction Reads in Spots Under Tissue","UMIs in Spots","Mean Reads per Spot","Median UMI Counts per Spot","Median Genes per Spot","Total Genes Detected","Reads Mapped to Genome","Reads Mapped Confidently to Genome","Reads Mapped Confidently to Intergenic Regions","Reads Mapped Confidently to Intronic Regions","Reads Mapped Confidently to Exonic Regions"]:
    key2 = trans.get(key,key)
    metrics_summary.append([key,str(ms[key2]).replace(',','')])
print(pd.DataFrame(metrics_summary))
pd.DataFrame(metrics_summary).T.to_csv(MetricsSummary,index=False,header=False)



f = open(Temp,encoding='utf-8').read()
for key in total_summary:
    for info in total_summary[key]:
        f = f.replace(f"$${info}$$",str(total_summary[key][info]))
        
with open(output,"w",encoding='utf-8') as out:
    out.write(f)
os.system(f"cp -r {SoftHome}/template/static {AnalysisDir}/{sample}/outs")
os.system(f"mkdir {AnalysisDir}/DynamicST&&ln -s {AnalysisDir}/{sample} {AnalysisDir}/DynamicST/")


