from dockstring import load_target
import csv
import argparse
import rdkit
from rdkit import Chem

def strip_mol(smi):
    """ if the smiles has a "." in there, usually meaning its a salt or
    multicomponent mol, return the biggest one (most characters). vina cant
    handle multidock so this is necessary
    probably replace this with the rdkit remove salt thing"""
    if smi.count(".")>0: #take heaviest mol if there are more
        bestsmi=""
        for csmi in smi.split("."):
            if len(csmi)>len(bestsmi):
                bestsmi=csmi
        smi=bestsmi
    return smi

def load_and_split(filenames, curr_slice, total_slices):
    """ load the input file(s), if it ends in sdf its treated as sdf
    otherwise its treated as a smi/csv of smiles. then, split the data
    (this is done so i can easily run a lot of jobs in parallel)"""
    cmpds_smis=[]
    for filename in filenames:
        if filename[-3:]=="sdf":
            suppl = Chem.SDMolSupplier(filename)
            for mol in suppl:
                cmpds_smis.append(strip_mol(Chem.MolToSmiles(mol)))
        else:
            with open(filename, newline='') as csvfile:
                reader = csv.reader(csvfile, delimiter=',', quotechar='"')
                for row in reader:
                    cmpds_smis.append(row[0])
    i1=int((curr_slice-1)/total_slices*len(cmpds_smis))
    i2=int(curr_slice/total_slices*len(cmpds_smis))
    return cmpds_smis[i1:i2]
    
def dock(cmpds,target,output):
    """use dockstring to dock the slice of compounds
    and then add Hs and write to SDF
    """
    mols = []
    scores = []
    for smi in cmpds:
        try:
            score,aux = target.dock(smi,pH=7.4)
            scores.append(score)
            mols.append(Chem.AddHs(aux["ligand"],addCoords=True))
        except:
            print("docking failure for smi {}".format(smi))
    writer = Chem.SDWriter(output)        
    for mol in mols:
        writer.write(mol)
    return scores



if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Fragment docking script.')
    parser.add_argument('--input', metavar='i',nargs='+', type=str, help='path to file(s) with to be docked fragments. can be sdf, can be smi')
    parser.add_argument('--target', metavar='t', required=False, type=str, default="ABL1", help='target to dock to (input must be compatible with dockstring)')
    parser.add_argument('--current_slice', metavar='c', required=False, type=int, default=1, help='if dataset is split, which slice to process. index count starts at 1')
    parser.add_argument('--slices', metavar='s', required=False, type=int, default=1, help='split dataset in n slices')
    parser.add_argument('--mode', metavar='m', required=False, type=str, default="dock", help='mode for generating conformers')    
    parser.add_argument('--output', metavar='o', type=str, required=False, default="output_docked.sdf", help='output sdf file with merged fragments')    
    args = parser.parse_args()
    target = load_target(args.target)
    cmpds = load_and_split(args.input,args.current_slice,args.slices)
    print("docking {} compounds".format(len(cmpds)))
    dock(cmpds,target,args.output)
   
