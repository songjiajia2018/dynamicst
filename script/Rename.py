import pandas as pd
import sys
import os
from pathlib import Path

starsummary,wpssummary,threads,pigz = sys.argv[1:5]

df = pd.read_csv(starsummary,index_col=0,header=None)
df.index = df.index.str.replace("Cells","Spots").str.replace("Cell","Spot")
df = df.dropna().T
df.to_csv(wpssummary,index=False)

filtered = Path(starsummary).parent/"filtered"
raw = Path(starsummary).parent/"raw"

if not list(filtered.glob("*gz")):
    os.system(f"{pigz} {threads} {filtered}/*")
    
if not list(raw.glob("*gz")):
    os.system(f"{pigz} {threads} {raw}/*")
