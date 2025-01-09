#python ../../smiles2sdf.py -i ligprep_pfkfb3.smi --align --heavy-atom-only true  --complete-rings-only true  --ring-compare StrictRingFusion --atom-compare CompareElements --bond-compare CompareOrderExact -o test_heavyatom_only

python ../../../smiles2sdf.py -i ligprep_pfkfb3.smi --align -r reference-ligand.sdf --heavy-atom-only true  --complete-rings-only true  --ring-compare StrictRingFusion --atom-compare CompareAnyHeavyAtom --bond-compare CompareOrderExact 

cat ligprep_pfkfb3_*sdf > combine.sdf 

