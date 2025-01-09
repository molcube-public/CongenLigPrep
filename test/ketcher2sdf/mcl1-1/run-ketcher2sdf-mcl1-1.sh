python ../../../ketcher2sdf.py -i ketcher-mcl1-1.mol --align --heavy-atom-only true  --complete-rings-only true  --ring-compare StrictRingFusion --atom-compare CompareElements --bond-compare CompareOrderExact  

cat ketcher-mcl1-1_*.sdf > combine.sdf 

