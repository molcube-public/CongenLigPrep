import os
import sys
import tempfile
import shutil  # for directory removal
import argparse
#sys.path.append('./libs')
sys.path.append(os.path.join(os.path.dirname(__file__), 'libs'))
from core_2_utils import *


def ketcher_combination_to_sdf(input_mdl_file, output_dir=None):
    """
    Convert Ketcher combination file to SDF files
    """ 
    if output_dir is None:
        output_dir = input_mdl_file[:-4]


    # create output directory if not exists 
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Create temp directory without context manager
    temp_dir = tempfile.mkdtemp()

    # parse MDL file into main structure and R-groups
    main_mol, r_groups = parse_mdl_file(input_mdl_file)

    # assemble molecules
    combined_ligs = assemble_molecule(main_mol, r_groups) 

    # write to SDF files
    output_sdf_files = []
    for i in range(len(combined_ligs)):

        write_mol_file(combined_ligs[i][0], f'{temp_dir}/ligand-{i+1}.sdf', combined_ligs[i][1])
        
        # one round of optimization
        #rdkit_optimize_molecule(f'{temp_dir}/ligand-{i+1}.sdf', f'{output_dir}/ligand-{i+1}.sdf')
          
        # two rounds of optimization
        rdkit_optimize_molecule(f'{temp_dir}/ligand-{i+1}.sdf', f'{temp_dir}/ligand-{i+1}_r1.sdf')
        rdkit_optimize_molecule(f'{temp_dir}/ligand-{i+1}_r1.sdf', f'{output_dir}/ligand-{i+1}.sdf')


        output_sdf_files.append(f'{output_dir}/ligand-{i+1}.sdf')
        if False:
            # check main structure atom indices in assembled molecule
            print(combined_ligs[i].main_str_atom_index) 

    # clean up temp directory
    shutil.rmtree(temp_dir)

    return output_sdf_files





if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert Ketcher combinatortial ligand MDL Molfile to SDF files')
    parser.add_argument('-i', '--input_mdl_file', 
                       required=True,
                       help='Input MDL Molfile (v2000) path')
    parser.add_argument('-o', '--output_dir', 
                       help='Output directory to save SDF files (default: same as input file without extension)',
                       default=None)
    
    args = parser.parse_args()
    output_dir = args.output_dir if args.output_dir else args.input_mdl_file[:-4]
    
    output_sdf_files = ketcher_combination_to_sdf(args.input_mdl_file, output_dir)
    #print(output_sdf_files)
