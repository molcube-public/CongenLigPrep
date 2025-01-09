 python ../../../smiles2sdf.py -i ligprep_syk.smi --align -r reference.sdf  --heavy-atom-only true --complete-rings-only true --atom-compare CompareAnyHeavyAtom  

cat ligprep_syk_*.sdf > combine.sdf

