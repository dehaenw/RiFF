import dock_frags
import rdkit
from rdkit import Chem
import csv
import argparse
from rdkit.Chem import RDConfig
import os
import sys
sys.path.append(os.path.join(RDConfig.RDContribDir, 'CalcLigRMSD'))
from CalcLigRMSD import *

def rmsd_docked(mol,dockedmol):
    """ function for calculating the rmsd between the atoms of the merged and
    docked poses
    """
    confmols=[]
    rmsds=[]
    for conf in dockedmol.GetConformers():
        molx=Chem.RWMol(dockedmol)
        molx.RemoveAllConformers()
        molx.AddConformer(conf)
        confmols.append(molx)
    for i,cm in enumerate(confmols):
        rmsd=CalcLigRMSD(cm,mol, rename_lig2 = False)
        rmsds.append(rmsd)
    return rmsds

def dock_compounds_and_generate_metrics(inputfile,target):
    """ open sdf file of merged compounds, redock and calculate rmsd between
    starting conformation and redocked conformation. also display the docking
    score
    """
    suppl = Chem.SDMolSupplier(inputfile)
    metrics = []
    for mol in suppl:
        try:
            smi = Chem.MolToSmiles(mol)
            print("docking", smi)
            score,aux = target.dock(smi,pH=7.4)
            metrics.append(([smi]+aux["affinities"]+rmsd_docked(mol,aux["ligand"])))
            print(metrics[-1])
        except Exception as e:
            try:
                print("docking failure for smi {}".format(smi))
                print(e)
            except:
                print("serious problem. couldnt convert mol in sdf to SMILES")
    return metrics
    
def save_report(metrics,outputfile):
    with open(outputfile, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',',quotechar='"', quoting=csv.QUOTE_MINIMAL)
        for metric in metrics:
            writer.writerow(metric)







if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Post merging docking validator.')
    parser.add_argument('--input', metavar='i',type=str, help='path to ONE file with the merged molecules, as sdf')
    parser.add_argument('--target', metavar='t', required=False, type=str, default="ABL1", help='target to dock to (input must be compatible with dockstring)')  
    parser.add_argument('--output', metavar='o', type=str, required=False, default="validation_report.csv", help='output with rmsd and docking score per compound')    
    args = parser.parse_args()
    target = dock_frags.load_target(args.target)
    metrics = dock_compounds_and_generate_metrics(args.input,target)
    save_report(metrics,args.output)

