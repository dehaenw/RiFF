
"""
clean this up
"""


import dock_frags

from rdkit import Geometry
from rdkit.Chem import ChemicalFeatures,rdDistGeom,rdMolTransforms
from rdkit.Chem.Pharm3D import Pharmacophore,EmbedLib
from rdkit.Numerics import rdAlignment
from rdkit import RDConfig
featFactory = AllChem.BuildFeatureFactory(RDConfig.RDDataDir+"/BaseFeatures.fdef")
from operator import itemgetter

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

feats1 = [ChemicalFeatures.FreeChemicalFeature('Aromatic', Geometry.Point3D(0.0, 0.0, 0.0)),
         ChemicalFeatures.FreeChemicalFeature('Donor', Geometry.Point3D(2.0, 0.0, 0.0)),
          ChemicalFeatures.FreeChemicalFeature('Acceptor', Geometry.Point3D(3.0, 1.0, 0.0)),]
ph4_1= Pharmacophore.Pharmacophore(feats1)
applyRadiiToBounds([1,1,1],ph4_1)
feats1 = [ChemicalFeatures.FreeChemicalFeature('Aromatic', Geometry.Point3D(0.0, 0.0, 0.0)),
         ChemicalFeatures.FreeChemicalFeature('Donor', Geometry.Point3D(2.0, 0.0, 0.0)),
          ChemicalFeatures.FreeChemicalFeature('Acceptor', Geometry.Point3D(3.0, 1.0, 0.0)),]
ph4_1= Pharmacophore.Pharmacophore(feats1)
applyRadiiToBounds([1,1,1],ph4_1)



mol = Chem.AddHs(Chem.MolFromSmiles("c1occc1OCCO"))
display(mol)

a,b = EmbedLib.MatchPharmacophoreToMol(mol,featFactory,ph4_1)
if a==True:
    boundsMat = rdDistGeom.GetMoleculeBoundsMatrix(mol)
    failed,boundsMatMatched,matched,matchDetails = EmbedLib.MatchPharmacophore(b,boundsMat,ph4_1)
    print(failed)
    if failed == 0:
        atomMatch = [list(x.GetAtomIds()) for x in matched]
        bm,embeddings,numFail = EmbedLib.EmbedPharmacophore(mol,atomMatch,ph4_1,count=100)
        if len(embeddings) > 0:
            SSD = TransformEmbeddings(ph4_1,embeddings,atomMatch)
            bestFitIndex = min(enumerate(SSD), key=itemgetter(1))[0] 
            ph4_score = SSD[bestFitIndex]
            conformation = embeddings[bestFitIndex]
        display(conformation)
        print(ph4_score)
        writer = Chem.SDWriter('tmp.sdf')
        writer.write(conformation)



if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Fragment joining script.')
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
   
