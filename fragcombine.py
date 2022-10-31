import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import argparse
import time

#note to be removed: why is this SMARTS so strange: because it needs to be able to account for both ammoniums and amines, as well as carboxylic acids and carboxylates. this is because during docking both states are possible for either.
reactionsdict = {"amide":{"r1":"[NX3H2-0,NX3H1-0,NX4H3+,NX4H2+;!$(NC=O);!$(NS=O);!$(NC=S):0]-[#1:1]","r2":"[CX3;H0:0](-[OX1H0-,OX2H1:1])(=[OX1])-[#6]","SMARTS":"[NX3H2-0,NX3H1-0,NX4H3+,NX4H2+;!$(NC=O);!$(NS=O);!$(NC=S):100](-[*:101]).[CX3;H0:102](-[OX1H0-,OX2H1:103])(=[OX1:104])-[#6:105]>>[NX3-0:100](-[*:101])-[*:102](=[*:104])-[*:105].[*:103]"},
"amine_alkylation":{"r1":"[NX3H2-0,NX3H1-0,NX4H3+,NX4H2+;!$(NC=O);!$(NS=O);!$(NC=S):0]-[#1:1]","r2":"[CX4:0]-[Cl,Br,I:1]","SMARTS":"[NX3H2-0,NX3H1-0,NX4H3+,NX4H2+;!$(NC=O);!$(NS=O);!$(NC=S):100](-[*:101]).[CX4:102]-[Cl,Br,I:103]>>[*:100](-[*:101])(-[*:102]).[*:103]"},
"suzuki":{"r1":"[#6:0]-[BX3;H0:1](-[OX2])(-[OX2])","r2":"[$([cX3](:*):*),$([cX2+](:*):*),$([CX3]=*),$([CX2+]=*):0]-[Cl,Br,I:1]","SMARTS":"[#6:100]-[BX3;H0:101](-[OX2])(-[OX2]).[$([cX3](:*):*),$([cX2+](:*):*),$([CX3]=*),$([CX2+]=*):102]-[Cl,Br,I:103]>>[*:100]-[*:102].[*:101].[*:103]"}}




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
    

    
def merge_mol(mol1,mol2,info1,info2,rxn):
    """ apply the reaction and make sure the mapped heavy atoms were in the
    mapping. explicit Hs are removed during the reaction (to make SMARTS
    more flexible) and added again afterwards.
    """
    merged_mols = rxn.RunReactants((Chem.RemoveHs(mol1),Chem.RemoveHs(mol2)))
    if len(merged_mols)>1:
        match_amounts=[]
        for curr_mol in merged_mols:
            product_idx = [int(atom.GetProp("react_atom_idx")) for atom in curr_mol[0].GetAtoms() if "old_mapno" in list(atom.GetPropNames())]
            match_amounts.append(product_idx.count(info1[0][0])+product_idx.count(info2[0][0]))
        merged_mol=merged_mols[np.argmax(match_amounts)][0] #get best match and account foridx collision
    else:
        merged_mol=merged_mols[0][0]
    Chem.SanitizeMol(merged_mol)
    return merged_mol
    
    
def distance_report(mol):
    """ make a report of:
    1. distances of atoms that are bonded
    2. distances of atoms with a topological distance 2-5 to eacht other
    3. distances of atoms with a topological distance 6+
    output is the minimum of category 1,2,3 and the average of the 5 lowest
    """
    noHs = mol
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
    
def combine_fragments(sdfs,reaction,outputpath,coord_cutoff=0.25,angle_cutoff=0.7,auto_threshold=False,auto_threshold_percentage=0.1,thr=[0.8,0.9,1.2,1.4,1.8,2.0]):
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
    rxn = AllChem.ReactionFromSmarts(reaction["SMARTS"])
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
    if auto_threshold==True:
        distances=[]
        angles=[]
        for i,i_r1c in enumerate(i_r1):
            for j,i_r2c in enumerate(i_r2):
                for i1 in i_r1c:
                    for i2 in i_r2c:
                        mse1,mse2 = calc_threshold(i1,i2)
                        distances.append(mse1)
                        angles.append(mse2)
        coord_argsort=np.argsort(distances)[:int(len(distances)*auto_threshold_percentage**0.5)]
        coord_cutoff=distances[coord_argsort[-1]]
        remaining_angles=[angles[x] for x in coord_argsort]
        angle_cutoff=np.sort(remaining_angles)[-int(len(remaining_angles)*auto_threshold_percentage**0.5)]
        print("autogenerated threshold at {} percent were coord {} and angle {}".format(100*auto_threshold_percentage,coord_cutoff,angle_cutoff))
    for i,i_r1c in enumerate(i_r1):
        for j,i_r2c in enumerate(i_r2):
            for i1 in i_r1c:
                for i2 in i_r2c:
                    mse1,mse2 = calc_threshold(i1,i2)
                    if mse1<coord_cutoff and mse2>angle_cutoff:
                        mmol = merge_mol(allmols[i],allmols[j],i1,i2,rxn)
                        if collision_detect(mmol,thr)==False:
                            winners.append(Chem.AddHs(mmol,addCoords=True))
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
    parser.add_argument('--input', metavar='i',nargs='+', type=str, help='path to sdf file(s) with to be joined fragment 3d conformations in there')
    parser.add_argument('--reaction', metavar='r', required=False, type=str, default="amide", help='reaction type from predefined dict')
    parser.add_argument('--custom_reaction', metavar='x',  required=False, type=str, help='custom reaction as reaction SMARTS. to be implemented')
    parser.add_argument('--coord_cutoff', metavar='c',  required=False, type=float, default=0.25, help='coord cutoff that determines overlap')
    parser.add_argument('--angle_cutoff', metavar='a',  required=False, type=float, default=0.7, help='angle cutoff that determines overlap')
    parser.add_argument('--auto_threshold', metavar='t',  required=False, type=bool, default=False, help='use automatic threshold instead of explicit one')
    parser.add_argument('--distance_cutoffs', metavar='d',nargs=6, required=False, type=float, default=[0.8,0.9,1.2,1.4,1.8,2.0], help='set distance based cutoffs') #try 1.0 1.2 2.0 2.2 2.5 3.0 for strict
    parser.add_argument('--auto_percentage', metavar='p',  required=False, type=float, default=0.01, help='percentage to use for thresolding. because of many post merge rejects it does not correspond t the final percentage')    
    parser.add_argument('--output', metavar='o', type=str, required=False, default="output.sdf", help='output sdf file with merged fragments')    
    args = parser.parse_args()
    combine_fragments(args.input,reactionsdict[args.reaction],args.output,args.coord_cutoff,args.angle_cutoff,args.auto_threshold,args.auto_percentage,args.distance_cutoffs)

