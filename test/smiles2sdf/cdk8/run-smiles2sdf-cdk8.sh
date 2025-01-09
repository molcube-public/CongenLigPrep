python ../../../smiles2sdf.py -i ligprep_cdk8.smi --align --heavy-atom-only true  --complete-rings-only true  --ring-compare StrictRingFusion --atom-compare CompareElements --bond-compare CompareOrderExact 

cat ligprep_cdk8_*.sdf > combine.sdf 


