import rdkit 
from rdkit import Chem
from rdkit.Chem import AllChem,Descriptors
import numpy as np
import csv
import os
import random

FRAGS_DIR = "path/to/frags"

"""
this is the script used to split and filter fragment databases
"""

amine_p = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O);!$(NS=O):0]")
cooh_p = Chem.MolFromSmarts("[CX3;H0:0](-[OX2;H1:1])(=[OX1])-[#6]")
halo_p = Chem.MolFromSmarts("[$([cX3](:*):*),$([cX2+](:*):*),$([CX3]=*),$([CX2+]=*):0]-[Cl,Br,I:1]")
boro_p = Chem.MolFromSmarts("[#6:0]-[BX3;H0:1](-[OX2])(-[OX2])")
sp3halo_p = Chem.MolFromSmarts("[CX4:0]-[Cl,Br,I:1]")

all_smis=[]
namedict={}
a_match=[]
b_match=[]
c_match=[]
h_match=[]
s_match=[]
frag_files= ["{}/{}".format(FRAGS_DIR,x) for x in os.listdir(FRAGS_DIR)]
for f in frag_files:
    with open(f, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=' ', quotechar='"')
        for row in reader:
            if len(all_smis)%10000==0:
                print(len(all_smis))
            smi=row[0]
            try:
                mol= Chem.MolFromSmiles(smi)
                if Descriptors.MolWt(mol)<250:
                    if len(list(mol.GetSubstructMatch(amine_p)))>0:
                        a_match.append(1)
                    else:
                        a_match.append(0)
                    if len(list(mol.GetSubstructMatch(boro_p)))>0:
                        b_match.append(1)
                    else:
                        b_match.append(0)
                    if len(list(mol.GetSubstructMatch(cooh_p)))>0:
                        c_match.append(1)
                    else:
                        c_match.append(0)
                    if len(list(mol.GetSubstructMatch(halo_p)))>0:
                        h_match.append(1)
                    else:
                        h_match.append(0)
                    if len(list(mol.GetSubstructMatch(sp3halo_p)))>0:
                        s_match.append(1)
                    else:
                        s_match.append(0)
                    smi = Chem.MolToSmiles(mol)
                    all_smis.append(smi)
                    namedict[smi]=row[1]
            except Exception as e:
                print(e)
                print("ERROR for round trip of ",smi)
                
a_smis=[]
b_smis=[]
c_smis=[]
h_smis=[]
s_smis=[]
ch_smis=[]
for i,smi in enumerate(all_smis):
    if a_match[i]+b_match[i]+c_match[i]+h_match[i]+s_match[i]==1:
        if a_match[i]==1:
            a_smis.append(smi)
        if b_match[i]==1:
            b_smis.append(smi)
        if c_match[i]==1:
            c_smis.append(smi)
        if h_match[i]==1:
            h_smis.append(smi)
        if s_match[i]==1:
            s_smis.append(smi)
    if c_match[i]+h_match[i]==2:
        ch_smis.append(smi)
a_smis=list(set(a_smis))
b_smis=list(set(b_smis))
c_smis=list(set(c_smis))
h_smis=list(set(h_smis))
s_smis=list(set(s_smis))
ch_smis=list(set(ch_smis))

with open('fragments/a_all.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter=' ',quotechar='"', quoting=csv.QUOTE_MINIMAL)
    for smi in a_smis:
        writer.writerow([smi,namedict[smi]])
        
with open('fragments/b_all.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter=' ',quotechar='"', quoting=csv.QUOTE_MINIMAL)
    for smi in b_smis:
        writer.writerow([smi,namedict[smi]])

with open('fragments/c_all.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter=' ',quotechar='"', quoting=csv.QUOTE_MINIMAL)
    for smi in c_smis:
        writer.writerow([smi,namedict[smi]])
        
with open('fragments/h_all.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter=' ',quotechar='"', quoting=csv.QUOTE_MINIMAL)
    for smi in h_smis:
        writer.writerow([smi,namedict[smi]])
        
with open('fragments/s_all.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter=' ',quotechar='"', quoting=csv.QUOTE_MINIMAL)
    for smi in s_smis:
        writer.writerow([smi,namedict[smi]])
        
with open('fragments/ch_all.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter=' ',quotechar='"', quoting=csv.QUOTE_MINIMAL)
    for smi in ch_smis:
        writer.writerow([smi,namedict[smi]])   
        
for n in [100,200,1000]:
    with open('fragments/a_{}.csv'.format(n), 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=' ',quotechar='"', quoting=csv.QUOTE_MINIMAL)
        for smi in random.sample(a_smis, n):
            writer.writerow([smi,namedict[smi]])
    with open('fragments/b_{}.csv'.format(n), 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=' ',quotechar='"', quoting=csv.QUOTE_MINIMAL)
        for smi in random.sample(b_smis, n):
            writer.writerow([smi,namedict[smi]])
    with open('fragments/c_{}.csv'.format(n), 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=' ',quotechar='"', quoting=csv.QUOTE_MINIMAL)
        for smi in random.sample(c_smis, n):
            writer.writerow([smi,namedict[smi]])
    with open('fragments/h_{}.csv'.format(n), 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=' ',quotechar='"', quoting=csv.QUOTE_MINIMAL)
        for smi in random.sample(h_smis, n):
            writer.writerow([smi,namedict[smi]])
    with open('fragments/s_{}.csv'.format(n), 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=' ',quotechar='"', quoting=csv.QUOTE_MINIMAL)
        for smi in random.sample(s_smis, n):
            writer.writerow([smi,namedict[smi]])
    with open('fragments/ch_{}.csv'.format(n), 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=' ',quotechar='"', quoting=csv.QUOTE_MINIMAL)
        for smi in random.sample(ch_smis, n):
            writer.writerow([smi,namedict[smi]])
                                                              
