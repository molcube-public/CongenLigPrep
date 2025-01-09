#python ../../smiles2sdf.py -i ligprep_p38.smi --align --heavy-atom-only true --atom-compare CompareAnyHeavyAtom --bond-compare CompareOrder  

python ../../../smiles2sdf.py -i ligprep_p38.smi --align --heavy-atom-only true  --atom-compare CompareAnyHeavyAtom 

 cat ligprep_p38_*.sdf > combine.sdf


