#!/home/pandunhuang/anaconda3/envs/mappy/bin/python
#from Bio import SeqIO
import sys
import pandas as pd
import mappy as mp
import argparse
import os


parser = argparse.ArgumentParser()
parser.add_argument('--r1',required=True,help="R1 fastq")
parser.add_argument('--r2',required=True,help="R2 fastq")
parser.add_argument('--bx',required=True,help="barcodex.txt")
parser.add_argument('--by',required=True,help="barcodey.txt")
parser.add_argument('--name',required=True,help="sample name")

parser.add_argument("--transcriptome",default="/disk/reference/WPSRanger_ref/mm10/star/",help="Path of folder containing transcriptome reference,default:/disk/reference/WPSRanger_ref/mm10/star/")

args = parser.parse_args()

R1              = args.r1
R2              = args.r2
barcodex        = args.bx
barcodey        = args.by
transcriptome   = args.transcriptome
name            = args.name

lbc = len([i for i in open(barcodex)])**2
linker_seq = "AGGCCAGAGCATTCGATCCACGTGCTTGAG"


def mismacth(adapter):
    bases = set("ATCGN")
    mismacth_set = set()
    mismacth_set.add(adapter)
    for index,base in enumerate(adapter):
        for elem in set(base)^bases:
            adapter2 = list(adapter)
            adapter2[index] = elem
            mismacth_set.add("".join(adapter2))
    return mismacth_set
    
linker_set = mismacth(linker_seq)
polyT_set  = mismacth("TTTTTT")
total_reads          = 0
correct_linker_count = 0
number_of_polyT      = 0 

d={}
for ID,seq,Q,Tag in mp.fastx_read(R1,read_comment=True):
    total_reads += 1
    bx = seq[:8]
    linker = seq[8:38]
    if linker in linker_set:
        correct_linker_count += 1
    by = seq[38:46]
    umi = seq[46:58]
    polyT = seq[58:58+6]
    if polyT in polyT_set:
        number_of_polyT += 1
    barcode = bx+by
    d.setdefault(barcode,[])
    d[barcode].append(umi)

dd = {}
for barcode in d:
    dd.setdefault(barcode,{})
    dd[barcode]['umi counts'] = len(d[barcode])
    dd[barcode]['umi 种类'] = len(set(d[barcode]))
    dd[barcode]['umi counts占总reads比例'] = f"{len(d[barcode])/total_reads:.2%}"

df = pd.DataFrame(dd).T
df = df.sort_values(by="umi 种类",ascending=False)
df3000 = df.head(lbc+1000)
#df.head(3000).to_csv("reads_summary.csv")


df2 = pd.DataFrame([total_reads,f"{correct_linker_count}({correct_linker_count/total_reads:.2%})",f"{number_of_polyT}({number_of_polyT/total_reads:.2%})"]).T
df2.columns= ['Total Reads',"Linker总数(允许1错配)","PloyT总数(允许1错配)"]
df2.to_excel("total_summary.xlsx")


shell = f"STAR  --genomeDir {transcriptome} --readFilesCommand pigz -dc --clipAdapterType CellRanger4 --readFilesIn {R2} {R1} --limitGenomeGenerateRAM 100000000000 --limitBAMsortRAM 100000000000 --soloType CB_UMI_Complex --soloCBwhitelist {barcodex} {barcodey} --soloBarcodeReadLength 0 --runThreadN 35 --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --soloOutFileNames outs --soloCellFilter TopCells 10000 --outSAMtype BAM SortedByCoordinate --outSAMattributes GX GN CB UB --soloCBmatchWLtype 1MM --soloFeatures Gene --outSAMunmapped Within --outSAMmultNmax 1 --soloCBposition 0_0_0_7 0_38_0_45 --soloUMIposition 0_46_0_57"

os.system(shell)
os.system("sambamba view -t 20 Aligned.sortedByCoord.out.bam -o Aligned.sortedByCoord.out.sam")

f = open("Aligned.sortedByCoord.out.sam")
d = {}

for i in f:
    i = i.strip().split("\t")
    CB,UB = i[-2:]
    CB = CB.replace("_","")[5:]
    d.setdefault(CB,[])
    d[CB].append(UB)

df_corr = pd.DataFrame(d.items())
df_corr['矫正后umi counts'] = df_corr[1].apply(lambda e:len(e))
df_corr['矫正后umi 种类']   = df_corr[1].apply(lambda e:len(set(e)))
df_corr['矫正后umi counts占总reads比例'] = df_corr['矫正后umi counts'].apply(lambda e:f"{e/total_reads:.2%}")
df_corr = df_corr.drop(1,axis=1)
df_corr = df_corr.set_index(0)
df_corr = df_corr[df_corr.index != "-"]

pd.concat([df3000,df_corr],axis=1).to_excel(f"{name}_矫正前后对比.xlsx")



from matplotlib import colors
import matplotlib.pyplot as plt
plt.rcParams["figure.facecolor"] = "white"

bcounts = df_corr.to_dict()['矫正后umi 种类']
def plot():
    cmap = colors.LinearSegmentedColormap.from_list("..",["#5658A6","#7BCAA4","#FDFEBD","#F99454","#CD2626"],N=1000)###自定义colorbar
    #df = pd.read_csv("/disk/pipeline/DynamicST/db/barcode_coordinate.txt",sep='\t',index_col=0)
    df = pd.DataFrame([[f"{bx.strip()}{by.strip()}",xindex,yindex] for xindex,bx in enumerate(open(barcodex,encoding='utf-8'),1) for yindex,by in enumerate(open(barcodey,encoding='utf-8'),1)])
    df.columns = ['barcode','x','y']
    df = df.set_index("barcode")
    xvalue = []
    yvalue = []
    color = []
    plt.figure(figsize=(16,14),facecolor="white",dpi=100)
    for barcode,(x,y) in df.iterrows():
        if barcode in bcounts:
            xvalue.append(x)
            yvalue.append(y)
            color.append(bcounts.get(barcode,0))
        
    plt.scatter(xvalue,yvalue,c=color,cmap =cmap,marker='s',s=80,linewidths=1,edgecolors="black")
    plt.ylim(max(df["x"])+1,0)
    plt.xlim(0,max(df["x"])+1)
    plt.axis("off")
    plt.colorbar(shrink=0.7)
    plt.subplots_adjust(bottom=0,top=1,left=0,right=1)  
    plt.savefig(f"{name}_umi_counts_distribution.png")
plot()
