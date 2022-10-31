# RiFF (working repo)
here's the riff working repo
## fragments
i have fragments categorized per functional group in the folder fragments
## fragcombine
the main script, heres an example usage
```
python3 fragcombine.py --i example_data/c_abl1_dock.sdf example_data/a_abl1_dock.sdf --angle_cutoff 0.85 --coord_cutoff 0.5 --output example_data/abl1_dock_merged.sdf
```
here's the result, 76/1000000 are hits, in the gif below the combined fragments are sshown in cyan:


![pic of the combined fragments inside the abl1 pocket](pictures/abl_riff.gif)
## dock frags
this is just a handy script for running dockstring jobs. the example data in example_data was generated using
```
python3 dock_frags.py --input fragments/c_1000.csv --target ABL1 --output example_data/c_abl1_dock.sdf
python3 dock_frags.py --input fragments/a_1000.csv --target ABL1 --output example_data/a_abl1_dock.sdf
```
## pharmacophore embed
provide a pharmacophore as an csv file of the form feature,x,y,z,radius, one line one feature.
for example:
```
Acceptor,2.275,0.0,0.0,1
Donor,-2.275,0.0,0.0,1
Aromatic,0.0,2.275,2,3.85,0.8
Acceptor,0.0,-2.275,3.85,1
```
embed fragment using
```
python3 pharmacophore_embed.py --input fragments/a_1000.csv --pharmacophore example_data/ph2.csv --retained_features 1 2 --output example_data/ph4_a_1000.sdf
python3 pharmacophore_embed.py --input fragments/a_1000.csv --pharmacophore example_data/ph2.csv --retained_features 3 4 --output example_data/ph4_a_1000.sdf
```
