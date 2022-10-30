import csv
import rdkit
from rdkit import Chem

def strip_mol(smi):
    """ if the smiles has a "." in there, usually meaning its a salt or
    multicomponent mol, return the biggest one (most characters). vina as well
    as other software cant handle multidock so this is necessary
    !!probably replace this with the rdkit remove salt thing"""
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
    
    
if __name__=="__main__":
    assert strip_mol("CCCCCN.HCl")=="CCCCCN"
