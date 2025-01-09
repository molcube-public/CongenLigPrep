#python ../../smiles2sdf.py -i ligprep_hif2a.smi --align --heavy-atom-only true  --complete-rings-only true  --ring-compare IgnoreRingFusion --atom-compare CompareElements --bond-compare CompareOrderExact
python ../../../smiles2sdf.py -i ligprep_hif2a.smi --align -r reference-ligand.sdf --heavy-atom-only true --match-chiral-tag false  --complete-rings-only false --ring-matches-ring-only false  --ring-compare IgnoreRingFusion --atom-compare CompareElements  --bond-compare CompareAny


cat ligprep_hif2a_*sdf  > combine.sdf 

