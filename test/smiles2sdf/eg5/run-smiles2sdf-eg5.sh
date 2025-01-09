python ../../../smiles2sdf.py -i ligprep_eg5.smi --align --heavy-atom-only true  --complete-rings-only true  --ring-compare StrictRingFusion --atom-compare CompareElements --bond-compare CompareOrderExact

 cat ligprep_eg5_*sdf > combine.sdf 

