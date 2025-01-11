# CongenLigPrep
![movie](https://github.com/user-attachments/assets/a0c9bb3c-ab4c-4021-a8cf-1ef89a7f93e9)




CongenLigPrep is an open-source utility to prepare high-quality congeneric series of ligand structures for alchemical free energy calculations.


## Feature and Components

![figures-1](https://github.com/user-attachments/assets/967e4a43-e49a-49de-bce6-62bba7a4df78)



CongenLigPrep consists of three Python scripts which accept different input file formats:

- sdf2sdf.py 
	- input: two ligand SDF files (target ligand and query ligand)
	- output: one ligand SDF file for aligned query ligand  
- ketcher2sdf.py
	- input file: MDL Molfile (V2000) for combinatorial ligands drawn on epam Ketcher
	- output: congeneric series of ligand SDF files
- smiles2sdf.py
	- input file: one SMILES file collecting multiple SMILES entries (containing H atoms explicitly)
	- output: congeneric series of ligand SDF files


## Installation

If you have miniconda installed, please run the following commands to install necessary dependencies and activate the new conda environment `congenligprep-1`.
```
conda env create -f environment.yml
conda activate congenligprep-1
```


## Quick Start

### Command to run 


sdf2sdf.py: 
```
python path/to/sdf2sdf.py -mol1 target.sdf -mol2 query.sdf -o query_aligned.sdf 
```

ketcher2sdf.py:
```
python path/to/ketcher2sdf.py -i ketcher.mol --align 
```

smiles2sdf.py: 
```
python path/to/smiles2sdf.py -i ligands.smi --align  
```

To see the complete option list, you may run
```
python path/to/any/script.py -h 
```


### Run all test examples

To run all test examples for ketcher2sdf.py: 
```
cd test/ketcher2sdf
bash run-all-test.sh
```

To run all test examples for smiles2sdf.py: 
```
cd test/smiles2sdf
bash run-all-test.sh
```



## Note about using epam Ketcher drawing 

epam Ketcher is an open-source software that allows you to draw 2D chemical structures in web browsers.

You may either use the online version (https://lifescience.opensource.epam.com/ketcher/demo.html) or the downloaded one (in ketcher-standalone-v2.26/main.html).

We prepared a tutorial video on how to draw combinatorial ligands with R-group tools for your reference. Please see https://www.youtube.com/watch?v=pretCk1KidI.




## Licensing 

CongenLigPrep is distributed under GPL-2.0 license. 
 

 
