python ../../../ketcher2sdf.py -i ketcher-mcl1-5.mol --align --heavy-atom-only true  --complete-rings-only true  --ring-compare StrictRingFusion --atom-compare CompareAnyHeavyAtom --bond-compare CompareOrderExact 

cat ketcher-mcl1-5_*.sdf > combined.sdf
