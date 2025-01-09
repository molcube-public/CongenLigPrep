python ../../../smiles2sdf.py -i ligprep_mcl1.smi --align -r 4hw3-native-ligand.sdf --heavy-atom-only true  --complete-rings-only true  --ring-compare IgnoreRingFusion --atom-compare CompareAnyHeavyAtom --bond-compare CompareOrderExact

cat ligprep_mcl1_*_aligned.sdf > combine_aligned.sdf
