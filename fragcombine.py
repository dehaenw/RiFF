import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import argparse
import time

reactionsdict = {"amide":{"r1":"[NX3,NX4+;H3,H2,H1;!$(NC=O);!$(NS=O):0]-[#1:1]","r2":"[CX3;H0:0](-[OX1H0-,OX2H1:1])(=[OX1])-[#6]","bondtype":Chem.BondType.SINGLE,"SMARTS":"[NX3,NX4+;H3,H2,H1;!$(NC=O);!$(NS=O):100]-[#1:101].[CX3;H0:102](-[OX1H0-,OX2H1:103])(=[OX1:104])-[#6:105]>>[*:100]-[*:102](=[*:104])-[*:105].[*:101]-[O:103]"}}



def get_coords_angles(mol,smartspattern):
    """ given a mol and SMARTS pattern, return the coords and angle (expressed
    as coord on the unit sphere) of the bond between atom with label 0 and 1.
    example:
    get_coords_angles(mol, '[NX3;H2,H1;!$(NC=O);!$(NS=O):0]-[#1:1]')
    for checking uncharged secondary and primary amines
    """
    match = mol.GetSubstructMatches(Chem.MolFromSmarts(smartspattern))
    matches = []
    for curr_match in match:
        conf = mol.GetConformer()
        #coord is the center of the bondwhere the merging will happen
        coords0 = tuple(conf.GetAtomPosition(curr_match[0]))
        coords1 = tuple(conf.GetAtomPosition(curr_match[1]))
        coords = tuple([0.5*coords0[i]+0.5*coords1[i] for i in range(3)])
        #direction as a point on the unit circle
        direction = [coords1[i]-coords0[i] for i in range(3)]
        a=sum([(x)**2 for x in direction])**0.5
        direction = tuple([x/a for x in direction])
        matches.append((curr_match[:2],coords,direction))
    return matches
    
def calc_threshold(info1,info2):
    """ given the info from coords_and_angles, calculate the deviation between
    two complementary groups."""
    coords_mse = sum([(info1[1][i]-info2[1][i])**2 for i in range(3)])/3
    direction_mse = sum([(info1[2][i]-info2[2][i])**2 for i in range(3)])/4 #in this case the range is [0-1]
    return coords_mse,direction_mse
    

    
def merge_mol(mol1,mol2,info1,info2,reactionsmarts):
    """ cut the bond definted by the initial tuple in info1 and info2 and 
    connect the first 2 atoms to make a merged molecule. needs the info from
    coords_and_angles.
    """
    rxn=AllChem.ReactionFromSmarts(reactionsmarts)
    merged_mols = rxn.RunReactants((mol1,mol2))
    merged_mol=None
    if len(merged_mols)>1:
        for curr_mol in merged_mols:
            old_idx = [int(atom.GetProp("react_atom_idx")) for atom in curr_mol[1].GetAtoms() if "old_mapno" in list(atom.GetPropNames())]
            if info1[0][1] == old_idx[0] and info2[0][1] == old_idx[1]:
                merged_mol=curr_mol[0]
    else:
        merged_mol=merged_mols[0][0]
    #delete this
    if merged_mol==None:
        print(Chem.MolToSmiles(mol1),Chem.MolToSmiles(mol2),info1,info2)
    return merged_mol
    
    
def distance_report(mol):
    """ make a report of:
    1. distances of atoms that are bonded
    2. distances of atoms with a topological distance 2-5 to eacht other
    3. distances of atoms with a topological distance 6+
    output is the minimum of category 1,2,3 and the average of the 5 lowest
    """
    noHs = Chem.RemoveHs(mol)
    distances_dict = {}
    dm = Chem.Get3DDistanceMatrix(noHs).flatten()
    tm = Chem.GetDistanceMatrix(noHs).flatten()
    i=1
    while np.count_nonzero(tm == i)>0:
        distances_dict[i]=[d for j,d in enumerate(dm) if tm[j]==i]
        i+=1
    bondmin = min(distances_dict[1])
    bondmin_mean5 = np.average(np.sort(distances_dict[1])[:5])
    closedistances = []
    i=2
    while i in distances_dict and i<=5:
        closedistances += distances_dict[i]
        i+=1
    closemin = min(closedistances)
    closemin_mean5 = np.average(np.sort(closedistances)[:5])
    fardistances = []
    i=6
    while i in distances_dict:
        fardistances += distances_dict[i]
        i+=1
    if len(fardistances)>0:
        farmin = np.min(fardistances)
        farmin_mean5 = np.average(np.sort(fardistances)[:5])
    else:
        farmin = None
        farmin_mean5 = None
    report = [bondmin,bondmin_mean5,closemin,closemin_mean5,farmin,farmin_mean5]
    return report

def collision_detect(mol,thr=[0.8,0.9,1.2,1.4,1.8,2.0]):
    """use the distance_report() output to detect if there is a collision or not"""
    dr = distance_report(mol)
    if dr[4]==None:
        dr[4:5]=[100,100]
    collision=True
    if dr[0]>thr[0] and dr[1]>thr[1] and dr[2]>thr[2] and dr[3]>thr[3] and dr[4]>thr[4] and dr[5]>thr[5]:
        collision=False
    return collision
    
def combine_fragments(sdfs,reaction,outputpath,coord_cutoff,angle_cutoff):
    """ load the mols in the sdfs, search for the complementary substructs,
    get their coords and angles, retain the ones mathcing within the defined
    threshold and then filter for collisions. save the remaining molecules as
    as an sdf file.
    """
    t1=time.time()
    winners = []
    i_r1 = []
    i_r2 = []
    r1_matches = 0
    r2_matches = 0
    allmols=[]
    #i do this be4hand and not in the loop bvecause otherwise i have to recalc the i_acids len(i_amines) times
    for curr_sdf in sdfs:
        suppl = Chem.SDMolSupplier(curr_sdf, removeHs=False)
        for mol in suppl:
            i_r1.append(get_coords_angles(mol, reaction["r1"]))
            if len(i_r1[-1])>0:
                r1_matches += 1
            i_r2.append(get_coords_angles(mol, reaction["r2"])) 
            if len(i_r2[-1])>0:
                r2_matches += 1
            allmols.append(mol)
    for i,i_r1c in enumerate(i_r1):
        for j,i_r2c in enumerate(i_r2):
            for i1 in i_r1c:
                for i2 in i_r2c:
                    mse1,mse2 = calc_threshold(i1,i2)
                    if mse1<coord_cutoff and mse2>angle_cutoff:
                        mmol = merge_mol(allmols[i],allmols[j],i1,i2,reaction["SMARTS"])
                        if collision_detect(mmol)==False:
                            winners.append(mmol)
    writer = Chem.SDWriter(outputpath)        
    for mmol in winners:
        writer.write(mmol)
    print("joining fragments took {} seconds".format(round(time.time()-t1,2)))
    print("total combined fragments is {} out of {} possibilities".format(len(winners),r1_matches*r2_matches))
    print("saved file to {}".format(outputpath))
    if len(winners)==0:
        print("i could not combine even one fragment. make sure the sdf has enough of each functional group that thresholds are ok and etc. make sure your sdfs have hydrogens in there.")
    return
    
    
if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Fragment joining script.')
    parser.add_argument('--input', metavar='-i',nargs='+', type=str, help='path to sdf file(s) with to be joined fragment 3d conformations in there')
    parser.add_argument('--reaction', metavar='-r', required=False, type=str, default="amide", help='reaction type from predefined dict')
    parser.add_argument('--custom_reaction', metavar='-x',  required=False, type=str, help='custom reaction as reaction SMARTS. to be implemented')
    parser.add_argument('--coord_cutoff', metavar='-c',  required=False, type=float, default=0.25, help='coord cutoff that determines overlap')
    parser.add_argument('--angle_cutoff', metavar='-a',  required=False, type=float, default=0.7, help='angle cutoff that determines overlap')
    parser.add_argument('--output', metavar='-o', type=str, required=False, default="output.sdf", help='output sdf file with merged fragments')    
    args = parser.parse_args()
    combine_fragments(args.input,reactionsdict[args.reaction],args.output,args.coord_cutoff,args.angle_cutoff)

