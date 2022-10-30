from dockstring import load_target
import csv
import argparse
import rdkit
from rdkit import Chem
import utils


def dock(cmpds,target,output,multiconf):
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
        if multiconf!=False:
            confs = mol.GetConformers()
            for c in confs:
                writer.write(mol, confId=c.GetId())
        else:
            writer.write(mol)
    return scores



if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Fragment docking script.')
    parser.add_argument('--input', metavar='i',nargs='+', type=str, help='path to file(s) with to be docked fragments. can be sdf, can be smi')
    parser.add_argument('--target', metavar='t', required=False, type=str, default="ABL1", help='target to dock to (input must be compatible with dockstring)')
    parser.add_argument('--current_slice', metavar='c', required=False, type=int, default=1, help='if dataset is split, which slice to process. index count starts at 1')
    parser.add_argument('--slices', metavar='s', required=False, type=int, default=1, help='split dataset in n slices')
    parser.add_argument('--multi_conformations', metavar='m', required=False, type=bool, default=False, help='so to true to retain all 9 dockstring conformatiosn.')    
    parser.add_argument('--output', metavar='o', type=str, required=False, default="output_docked.sdf", help='output sdf file with merged fragments')    
    args = parser.parse_args()
    target = load_target(args.target)
    cmpds = utils.load_and_split(args.input,args.current_slice,args.slices)
    print("docking {} compounds".format(len(cmpds)))
    dock(cmpds,target,args.output,args.multi_conformations)
   
