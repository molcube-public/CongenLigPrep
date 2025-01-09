python ../../../smiles2sdf.py  -i ligprep_thrombin.smi --align --match-chiral-tag false

cat ligprep_thrombin_*sdf >combine.sdf
