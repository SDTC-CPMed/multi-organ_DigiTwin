

import pandas as pd
import numpy as np
import scanpy as sc
from numpy.random import choice
from dca.api import dca
from scipy.stats import spearmanr
from sklearn.isotonic import IsotonicRegression
from dca.io import normalize
from collections import Counter
from pathlib import Path
import os
import sys

SSpath=sys.argv[1]
Refpath=sys.argv[2]
resolution=float(sys.argv[3])

ref =pd.read_csv(Refpath, index_col=0)
ref=ref.drop("FineCellType",1)
sc._settings.ScanpyConfig.figdir=Path("./data/clustering_and_celltype_identification")



os.makedirs("./data/clustering_and_celltype_identification", exist_ok=True)


RefLabels=ref.index.to_list()



ssData=pd.read_csv(SSpath, index_col=0)


ssData1=ssData.transpose()

#Merging data now

fullData=pd.concat([ssData1,ref], join='inner', sort=False)



nR=ref.shape[0]
nS=ssData1.shape[0]
res=sc.pp.filter_genes(fullData.iloc[:nS,:], min_counts=1, inplace=False)

p=fullData.shape[1]

fullData=fullData.iloc[:,res[0].to_list()].copy()

ind=choice(list(range(nS)), nS, replace=False)
ssData2=sc.AnnData(fullData.to_numpy()[ind,:])


model=dca(ssData2, verbose=True, return_model=True, return_info=True, hidden_dropout=0.01, epochs=300)


enc=model.get_encoder()
fullData2=fullData.copy()

fullData2.iloc[ind,:]=ssData2.X



def corr_spear(A,B):
    rho1,pi=spearmanr(A, B, axis=1)
    sA=A.shape[0]
    return rho1[:sA,sA:]



ssData3=fullData2.iloc[:nS,:].copy()
ref3=fullData2.iloc[nS:,:].copy()

rho0=corr_spear(np.log(ssData3.to_numpy()), ref3.to_numpy())

ssData4=fullData2.iloc[:nS,:].copy()
ref4=fullData.iloc[nS:,:].copy()

ref4n=ref4.copy()

nR=ref4.shape[0]


for i in range(nR):
    Cand=rho0[:,i]
    index=((-Cand).argsort())[:1]
    ssM=ssData4.iloc[index,:].mean(axis=0)

    
    x=ref4.iloc[i,:]
    y=np.log(ssM)
    iso_reg = IsotonicRegression().fit(x, y)
    yhat=iso_reg.predict(x)

    
    ref4n.iloc[i,:]=np.exp(yhat)



refData2=np.round(ref4n.to_numpy())
ind2=np.where(refData2.sum(axis=1)>0)

scData2=fullData.iloc[:nS,:].to_numpy()

full=sc.AnnData(np.concatenate((scData2, refData2)),  )

full1=normalize(full, filter_min_counts=False)

Xpred1=enc.predict({'count': full1.X, 'size_factors':full1.obs.size_factors })

ad=sc.AnnData(Xpred1)
ad.obs['Group']=['1']*nS+np.array(RefLabels)[ind2[0]].tolist()
sc.pp.neighbors(ad, n_neighbors=30)


sc.tl.leiden(ad, key_added="leiden_scvi", resolution=resolution)


Refclust=ad.obs["leiden_scvi"][nS:]
Sclust=ad.obs["leiden_scvi"][:nS]

Clusts=ad.obs["leiden_scvi"].unique()
HH={}

CellT={}
CellP={}

r4i=ref4.index[ind2[0]]
rhoS=rho0[:, ind2[0]]

for cl in Clusts:
    obs=np.where(Sclust==cl)[0]

    ii=np.where(Refclust==cl)[0]


    rfs=r4i[ii]

    wh=rfs.unique()

    k=0
    rhos=np.empty(shape=(len(wh),len(obs)))
    for i in wh:
        iii=np.where(rfs==i)[0]

        rhos[k,:]=rhoS[obs[:,None],ii[iii][None,:]].max(axis=1)#
        k=k+1


    if len(wh)>0 and len(obs)>0:

        CT=wh[list(Counter(rhos.argmax(axis=0)).keys())]
        Ws=np.array(list(Counter(rhos.argmax(axis=0)).values()))
        Prs=Ws/Ws.sum()
        Sumr=pd.DataFrame(list(zip(CT,Ws,Prs)),columns=["CellType", "Count", "Prob" ])
        CellT[cl]=CT[Prs.argmax()]
        CellP[cl]=Prs.max()

    else:
        Sumr=0
        CellT[cl]="Unknown"+str(cl)
        CellP[cl]=1

    HH[cl]=Sumr



Refclust=ad.obs["leiden_scvi"][nS:]
Sclust=ad.obs["leiden_scvi"][:nS]

Clusts=ad.obs["leiden_scvi"].unique()



Tlab=ad.obs["leiden_scvi"]
Labs=Tlab.to_list()
for i in range(len(Tlab)):
    Labs[i]=CellT[Tlab[i]]



#Processing to desired format
labsP=np.vstack(np.char.split(ssData1.index.to_numpy().astype(str), "_"))
Mouse=labsP[:,2].copy()
hi=np.where(labsP[:,2].astype('<U1')=="H")
si=np.where(labsP[:,2].astype('<U1')=="S")
labsP[hi,2]="Healthy"
labsP[si,2]="Sick"

labs=pd.DataFrame({"cell_ID":ssData1.index.to_numpy().astype(str),"CellType":Labs[:nS], 
                   "disease_group":labsP[:,2], "tissue": labsP[:,0], "sample":Mouse})


ad2=sc.AnnData(Xpred1[:nS,:])
sc.pp.neighbors(ad2, n_neighbors=30)
sc.tl.umap(ad2, min_dist=0.1)

ad2.obs['CellType']=labs["CellType"].tolist()


sc.pl.umap(ad2, color=["CellType"], save="clusters", show=False)


ad2.obs['Tissue']=labs["tissue"].tolist()
sc.pl.umap(ad2, color=["Tissue"], save="tissues", show=False)


ad2.obs['sample']=labs["sample"].tolist()
sc.pl.umap(ad2, color="sample", save="sample.pdf", show=False)


labs.to_csv('./data/clustering_and_celltype_identification/cluster_ids.csv')



