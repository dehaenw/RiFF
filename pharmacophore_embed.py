import utils
import rdkit
import argparse
import csv
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import Geometry
from rdkit.Chem import ChemicalFeatures,rdDistGeom,rdMolTransforms
from rdkit.Chem.Pharm3D import Pharmacophore,EmbedLib
from rdkit.Numerics import rdAlignment
from rdkit import RDConfig
featFactory = AllChem.BuildFeatureFactory(RDConfig.RDDataDir+"/BaseFeatures.fdef")
from operator import itemgetter
from rdkit import RDLogger                                                                                                                                                               
RDLogger.DisableLog('rdApp.*') 


import traceback
import sys



"""
note this code heavily makes use of 
https://github.com/rdkit/UGM_2016/blob/master/Notebooks/Stiefl_RDKitPh4FullPublication.ipynb
thanks to the authors.
"""
def applyRadiiToBounds(radii,pcophore):
    for i in range(len(radii)):
        for j in range(i+1,len(radii)):
            sumRadii = radii[i]+radii[j]
            pcophore.setLowerBound(i,j,max(pcophore.getLowerBound(i,j)-sumRadii,0))
            pcophore.setUpperBound(i,j,pcophore.getUpperBound(i,j)+sumRadii)
            
def GetTransformMatrix(alignRef,confEmbed,atomMatch):
    alignProbe = []
    for matchIds in atomMatch:
        dummyPoint = Geometry.Point3D(0.0,0.0,0.0)
        for id in matchIds:
            dummyPoint += confEmbed.GetAtomPosition(id)
        dummyPoint /= len(matchIds)
        alignProbe.append(dummyPoint)
    return (rdAlignment.GetAlignmentTransform(alignRef,alignProbe))

def TransformEmbeddings(pcophore,embeddings,atomMatch):
    alignRef = [f.GetPos() for f in pcophore.getFeatures()]
    SSDs = []
    for embedding in embeddings:
        conf = embedding.GetConformer()
        SSD,transformMatrix = GetTransformMatrix(alignRef,conf,atomMatch)
        rdMolTransforms.TransformConformer(conf,transformMatrix)
        SSDs.append(SSD)
    return(SSDs)
    
def try_embed(mol,ph4):
    conformation=None
    a,b = EmbedLib.MatchPharmacophoreToMol(mol,featFactory,ph4)
    if a==True:
        boundsMat = rdDistGeom.GetMoleculeBoundsMatrix(mol)
        failed,boundsMatMatched,matched,matchDetails = EmbedLib.MatchPharmacophore(b,boundsMat,ph4)
        if failed == 0:
            atomMatch = [list(x.GetAtomIds()) for x in matched]
            bm,embeddings,numFail = EmbedLib.EmbedPharmacophore(mol,atomMatch,ph4,count=100)
            if len(embeddings) > 0:
                SSD = TransformEmbeddings(ph4,embeddings,atomMatch)
                best10 = np.argsort(SSD)[:10]
                #bestFitIndex = min(enumerate(SSD), key=itemgetter(1))[0] 
                ph4_scores = [SSD[x] for x in best10]
                conformation = [embeddings[x] for x in best10]
                if ph4_scores[0]>1.0:
                    conformation=None
            else:
                conformation=None
    return conformation
    
def load_pharmacophore(path,retained_features=[-1]):
    """feed a csv with feature name,x,y,z,radius per line for every feature
    use retained features if you want to exclude some features
    """
    feature=1
    feats=[]
    radii=[]
    with open(path, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='"')
        for row in reader:
            if feature in retained_features or retained_features[0]==-1:
                feats.append(ChemicalFeatures.FreeChemicalFeature(row[0], Geometry.Point3D(float(row[1]), float(row[2]), float(row[3]))))
                radii.append(float(row[4]))
            feature += 1
    ph4= Pharmacophore.Pharmacophore(feats)
    applyRadiiToBounds(radii,ph4)                
    return ph4
    
def embed(cmpds,ph4,output):
    mols = []
    scores = []
    for smi in cmpds:
        try:
            mol = Chem.AddHs(Chem.MolFromSmiles(smi),addCoords=True)
            Chem.AssignStereochemistry(mol,cleanIt=True,force=True,flagPossibleStereoCenters=True)
            confmol = try_embed(mol,ph4)
            mols.append(confmol)
        except Exception as e:
            print("cant embed")
            print(e)
            print(traceback.format_exc())
            # or
            print(sys.exc_info()[2])
    writer = Chem.SDWriter(output)        
    for mol in mols:
        if mol!=None:
            for c in mol:
                writer.write(c)
    return


if __name__=="__main__":
    parser = argparse.ArgumentParser(description='pharmacophore embedding script.')
    parser.add_argument('--input', metavar='i',nargs='+', type=str, help='path to file(s) with to be embedded fragments. can be sdf, can be smi')
    parser.add_argument('--pharmacophore', metavar='p', required=False, type=str, default="example_data/ph1.csv", help='link to csv with feature+coords per line')
    parser.add_argument('--retained_features', metavar='r',nargs='+', type=int, default=[-1], help='which features are used (1-indexed). put to -1 to use all features')
    parser.add_argument('--current_slice', metavar='c', required=False, type=int, default=1, help='if dataset is split, which slice to process. index count starts at 1')
    parser.add_argument('--slices', metavar='s', required=False, type=int, default=1, help='split dataset in n slices')
    parser.add_argument('--output', metavar='o', type=str, required=False, default="output_embedded.sdf", help='output sdf file with merged fragments')    
    args = parser.parse_args()
    ph4 = load_pharmacophore(args.pharmacophore, args.retained_features)
    cmpds = utils.load_and_split(args.input,args.current_slice,args.slices)
    print("embedding {} compounds".format(len(cmpds)))
    embed(cmpds,ph4,args.output)
   
