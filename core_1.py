import rdkit
from rdkit.Chem import AllChem

import sys
import os
#sys.path.append('./libs')
sys.path.append(os.path.join(os.path.dirname(__file__), 'libs'))
from core_1_utils import *
from core_2_utils import rdkit_optimize_molecule
from gcutil import *

import tempfile
import shutil

import argparse

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

# this part is to read 2 ligand SDF files, find their MCS 



def align_with_zmat(mol1_sdf, mol2_sdf, processed_mol2_sdf, 
                         sym_equivalent_patterns=None, debug=False,
                         timeout=10, ring_matches_ring_only=True,
                         complete_rings_only=False, ring_compare='IgnoreRingFusion',
                         atom_compare='CompareElements', bond_compare='CompareOrderExact',
                         heavy_atom_only=False, optimize_h_atoms=False, zmat_gen_method=1, match_chiral_tag=True):
    """
    Process a pair of molecules to find MCS, identify substituents, and generate new coordinates.
    
    Args:
        mol1_sdf (str): Path to first molecule SDF file
        mol2_sdf (str): Path to second molecule SDF file
        processed_mol2_sdf (str, optional): Path to save processed second molecule
        sym_equivalent_patterns (list, optional): List of patterns for symmetry equivalent groups
            e.g. [{'anchor': 'C', 'neighbors': ['H', 'H', 'H']}, ...]
        debug (bool): Whether to print debug information
        timeout (int): Maximum time in seconds for MCS search
        ring_matches_ring_only (bool): If True, ring atoms only match ring atoms
        complete_rings_only (bool): If True, only complete rings are considered
        ring_compare (str): Ring matching behavior ('StrictRingFusion', 'IgnoreRingFusion', 'PermissiveRingFusion')
        atom_compare (str): Atom comparison method ('CompareAny', 'CompareElements', 'CompareIsotopes', 'CompareAnyHeavyAtom')
        bond_compare (str): Bond comparison method ('CompareAny', 'CompareOrder', 'CompareOrderExact')
        heavy_atom_only (bool): If True, only heavy atoms are considered during MCS search
        optimize_h_atoms (bool): If True, optimize H atoms; only available when heavy_atom_only is True
        zmat_gen_method (int): Method to generate z-matrix (1: use our own zmat generation, 2: use geomConvert's method, 3: use XYZ-to-ZMAT's method)
        match_chiral_tag (bool): If True, chiral tags are matched
    Returns:
        dict: Results containing processed molecules and coordinates
    """

    visualize_flag = False


    # Set default patterns if none provided
    if sym_equivalent_patterns is None:
        sym_equivalent_patterns = [
            {'anchor': 'C', 'neighbors': ['H', 'H', 'H']},  # -CH3
            {'anchor': 'N', 'neighbors': ['H', 'H']},      # -NH2
            {'anchor': 'N', 'neighbors': ['H', 'H', 'H']}, # -NH3^+
            {'anchor': 'C', 'neighbors': ['F', 'F', 'F']},  # -CF3
            {'anchor': 'C', 'neighbors': ['H', 'H', '?']},  # -CH2?
        ]

    # this part is to read 2 ligand SDF files, find their MCS 
    # 
    mol1 = rdkit.Chem.SDMolSupplier(mol1_sdf,removeHs=False,sanitize=True)[0]   
    mol2 = rdkit.Chem.SDMolSupplier(mol2_sdf,removeHs=False,sanitize=True)[0]   

    mol1_2d = rdkit.Chem.SDMolSupplier(mol1_sdf,removeHs=False,sanitize=True)[0]   
    mol2_2d = rdkit.Chem.SDMolSupplier(mol2_sdf,removeHs=False,sanitize=True)[0]  

    AllChem.Compute2DCoords(mol1_2d)
    AllChem.Compute2DCoords(mol2_2d)



    if heavy_atom_only:
        temp_dir = tempfile.mkdtemp()
        mol1_noH_temp_sdf = os.path.join(temp_dir, 'mol1_noH.sdf')
        mol2_noH_temp_sdf = os.path.join(temp_dir, 'mol2_noH.sdf')
        mol1_noH_sdf_text, mol1_noH_atom_mapping = remove_hydrogens_and_get_mapping_for_sdf(mol1_sdf)
        mol2_noH_sdf_text, mol2_noH_atom_mapping = remove_hydrogens_and_get_mapping_for_sdf(mol2_sdf)
        write_processed_sdf(mol1_noH_sdf_text, mol1_noH_temp_sdf)
        write_processed_sdf(mol2_noH_sdf_text, mol2_noH_temp_sdf)

        mol1_noH = rdkit.Chem.SDMolSupplier(mol1_noH_temp_sdf,removeHs=False,sanitize=True)[0]
        mol2_noH = rdkit.Chem.SDMolSupplier(mol2_noH_temp_sdf,removeHs=False,sanitize=True)[0]

        hits, r = find_mcs_matches_rdkit_2d(
            mol1_noH, mol2_noH,
            timeout=timeout,
            ring_matches_ring_only=ring_matches_ring_only,
            complete_rings_only=complete_rings_only,
            ring_compare=ring_compare,
            atom_compare=atom_compare,
            bond_compare=bond_compare,
            match_chiral_tag=match_chiral_tag,                                 
                                        )
        
        # use the mapping info in mol1_noH_atom_mapping/mol2_noH_atom_mapping to map hits to mol1/mol2
        hits_orig_1 = tuple(mol1_noH_atom_mapping[idx] for idx in hits[0])
        hits_orig_2 = tuple(mol2_noH_atom_mapping[idx] for idx in hits[1])
        hits = [hits_orig_1, hits_orig_2]
        shutil.rmtree(temp_dir)


    else: # heavy_atom_only is False
    # find MCS here 
        hits, r = find_mcs_matches_rdkit_2d(
            mol1, mol2,
            timeout=timeout,
            ring_matches_ring_only=ring_matches_ring_only,
            complete_rings_only=complete_rings_only,
            ring_compare=ring_compare,
            atom_compare=atom_compare,
            bond_compare=bond_compare,                                 
            match_chiral_tag=match_chiral_tag )



    if debug>1:
        print(hits[0])
        print(hits[1])
        display(rdkit.Chem.Draw.MolsToGridImage([mol1_2d, mol2_2d], highlightAtomLists=hits,subImgSize=(500,500) ) )

    # find substituents using adjacent matrix and connected_components 

    # For first ligand
    substituent_indices_1 = find_substituent_indices(mol1, hits[0])

    # For second ligand
    substituent_indices_2 = find_substituent_indices(mol2, hits[1])

    # Print results if needed
    if debug:
        for mol_num, indices in enumerate([substituent_indices_1, substituent_indices_2], 1):
            print(f"\nSubstituents for molecule {mol_num}:")
            print('Number of substituent sites:',len(indices))
            for site_label, atoms in indices.items():
                print(f"Substituent site {site_label} contains atoms: {atoms}")





    # Find anchor and substituent atoms for each substituent
    results_mol1 = collect_substituent_info(mol1, mol1_2d, substituent_indices_1, hits[0], visualize=visualize_flag, debug=debug)
    results_mol2 = collect_substituent_info(mol2, mol2_2d, substituent_indices_2, hits[1], visualize=visualize_flag, debug=debug)

    


    # Refine MCS using pattern information
    old_hits = list(hits)
    while True and not heavy_atom_only:
        # update hits for mol1
        hits = exclude_substituents_from_mcs(mol1, results_mol1, hits, sym_equivalent_patterns, mol_index=0)

        # update hits for mol2
        hits = exclude_substituents_from_mcs(mol2, results_mol2, hits, sym_equivalent_patterns, mol_index=1)


        # keep only the largest MCS 
        am_1 = Chem.GetAdjacencyMatrix(mol1)
        am_1_mcs = am_1[np.ix_(hits[0],hits[0])] 
        current_mcs_idx_largest = find_most_frequent_elements(connected_components(am_1_mcs)[1])
        new_hit_1 = get_sublist(hits[0], current_mcs_idx_largest)
        new_hit_2 = get_sublist(hits[1], current_mcs_idx_largest)
        hits = [new_hit_1, new_hit_2]

        # compare hits with old_hits, if no change, break the loop 
        if  all(sorted(h) == sorted(oh) for h, oh in zip(hits, old_hits)):    
            break
        else:
            old_hits = list(hits)


        if debug>1:
            display(rdkit.Chem.Draw.MolsToGridImage([mol1_2d, mol2_2d], highlightAtomLists=hits,subImgSize=(500,500) ) )



        # after MCS refinement 
        # find substituents using adjacent matrix and connected_components 

        # for first ligand
        substituent_indices_1 = find_substituent_indices(mol1, hits[0])

        # for second ligand
        substituent_indices_2 = find_substituent_indices(mol2, hits[1])

        # Print results if needed
        if debug:
            for mol_num, indices in enumerate([substituent_indices_1, substituent_indices_2], 1):
                print(f"\nSubstituents for molecule {mol_num}:")
                print('Number of substituent sites:',len(indices))
                for site_label, atoms in indices.items():
                    print(f"Substituent site {site_label} contains atoms: {atoms}")


        # Find anchor and substituent atoms for each substituent
        results_mol1 = collect_substituent_info(mol1, mol1_2d, substituent_indices_1, hits[0], visualize=visualize_flag, debug=debug)
        results_mol2 = collect_substituent_info(mol2, mol2_2d, substituent_indices_2, hits[1], visualize=visualize_flag, debug=debug)

    # generate z-matrix for substituents for mol2 

    m2sub_xyz_zmat_list = [] 
    for site_id, site_data in results_mol2.items():
        m2sub_xyz_zmat = extract_anchor_subst_xyz(mol2, site_data) 

        if zmat_gen_method == 1 or zmat_gen_method == 2:
            generate_zmatrix_file(m2sub_xyz_zmat, rvar=False, avar=False, dvar=False, zmat_gen_method=1)
            
        elif zmat_gen_method == 3:
            generate_zmatrix_file_engine2(m2sub_xyz_zmat)
        else:
            raise ValueError('Invalid zmat generation method')

        if debug>1: 
            print('zmat file:', m2sub_xyz_zmat['zmat_file'])
            print(m2sub_xyz_zmat['indices'][:3])       

        m2sub_xyz_zmat_list.append(m2sub_xyz_zmat)





    # build substituent XYZ coordinates for mol2 using z-matrix
    for m2sub_xyz_zmat in m2sub_xyz_zmat_list:

        m2sub_xyz_zmat = process_substituent_coord_xyz_zmat(m2sub_xyz_zmat, mol1, hits) 



    # Final update of coordinates for mol2
    # get XYZ coordinates of MCS atoms in mol1
    mol1_xyz = read_sdf_coordinates(mol1_sdf)

    mol1_mcs_xyz = mol1_xyz[list(hits[0])]


    orig_coords = read_sdf_coordinates(mol2_sdf)
    new_coords = orig_coords.copy()

    # replace MCS atoms
    new_coords = update_coordinates(new_coords, hits[1], mol1_mcs_xyz)

    # For each substituent
    for m2sub_xyz_zmat in m2sub_xyz_zmat_list:
        # Get indices and new coordinates from m2sub_xyz_zmat
        indices = m2sub_xyz_zmat['indices']
        subst_coords = m2sub_xyz_zmat['coordinates_new']  # or whatever key contains the new coordinates
        
        # Update coordinates
        new_coords = update_coordinates(new_coords, indices, subst_coords)


    # Export processed SDF file, optimize H atoms if needed
    if heavy_atom_only and optimize_h_atoms:
        temp_dir = tempfile.mkdtemp()
        mol2_temp_sdf = os.path.join(temp_dir, 'mol2_before_h_opt.sdf')
        new_sdf_file_tmp = modify_sdf_coordinates(mol2_sdf, new_coords, mol2_temp_sdf)
        rdkit_optimize_molecule(mol2_temp_sdf, processed_mol2_sdf, optimize_h_only=True)
        shutil.rmtree(temp_dir)

    else:
        new_sdf_file = modify_sdf_coordinates(mol2_sdf, new_coords, processed_mol2_sdf)

    return 



if __name__ == '__main__':
    import json

    class ArgumentParserWithValidation(argparse.ArgumentParser):
        def parse_args(self):
            args = super().parse_args()
            if args.optimize_h_atoms and not args.heavy_atom_only:
                self.error("--optimize-h-atoms can only be used when --heavy-atom-only is True")
            return args

    parser = ArgumentParserWithValidation(description='Align molecule B to molecule A using Z-matrix method')
    
    # Required arguments with named parameters
    parser.add_argument('-mol1', '--molecule1', required=True,
                        help='Path to target molecule SDF file')
    parser.add_argument('-mol2', '--molecule2', required=True,
                        help='Path to query molecule SDF file')
    
    # Optional arguments
    parser.add_argument('--output', '-o', 
                        help='Output SDF file path for processed molecule (default: [mol2]_processed.sdf)')
    parser.add_argument('--patterns', '-p',
                        help='Path to JSON file containing symmetry equivalent patterns')
    parser.add_argument('--timeout', '-t', type=int, default=10,
                        help='Maximum time in seconds for MCS search (default: 10)')
    parser.add_argument('--debug', '-d', type=int, default=0, 
                        help='Debug level (0: disabled, 1: basic info, 2: basic + 2D molecule visualization) (default: 0)')
    parser.add_argument('--ring-matches-ring-only', type=str2bool, default=True, metavar='BOOL',
                        help='Ring atoms only match ring atoms (default: True)')
    parser.add_argument('--complete-rings-only', type=str2bool, default=False, metavar='BOOL',
                        help='Only complete rings are considered (default: False)')
    parser.add_argument('--ring-compare', choices=['StrictRingFusion', 'IgnoreRingFusion', 'PermissiveRingFusion'],
                        default='IgnoreRingFusion',
                        help='Ring matching behavior (default: IgnoreRingFusion)')
    parser.add_argument('--atom-compare', choices=['CompareAny', 'CompareElements', 'CompareIsotopes', 'CompareAnyHeavyAtom'],
                        default='CompareElements',
                        help='Atom comparison method (default: CompareElements)')
    parser.add_argument('--bond-compare', choices=['CompareAny', 'CompareOrder', 'CompareOrderExact'],
                        default='CompareOrder',
                        help='Bond comparison method (default: CompareOrder)')
    parser.add_argument('--heavy-atom-only', type=str2bool, default=False, metavar='BOOL',
                        help='Only heavy atoms are considered during MCS search (default: False)')
    parser.add_argument('--optimize-h-atoms', type=str2bool, default=False, metavar='BOOL',
                        help='Optimize H atoms (default: False)')

    
    parser.add_argument('--zmat-gen-method', 
                      choices=[1, 2, 3], 
                      type=int,
                      default=1,
                      help='Method to generate z-matrix (1: use our own zmat generation, 2: use geomConvert\'s method, 3: use XYZ-to-ZMAT\'s method) (default: 1)')
    
    parser.add_argument('--match-chiral-tag', type=str2bool, default=True, metavar='BOOL',
                        help='Match chiral tags (default: True)')

    args = parser.parse_args()

    # Set default output path if not provided
    if args.output is None:
        args.output = args.molecule2[:-4] + '_processed.sdf'

    # Load patterns from file if provided, otherwise use defaults
    if args.patterns:
        with open(args.patterns, 'r') as f:
            sym_equivalent_patterns = json.load(f)
    else:
        sym_equivalent_patterns = [
            {'anchor': 'C', 'neighbors': ['H', 'H', 'H']},  # -CH3
            {'anchor': 'N', 'neighbors': ['H', 'H']},      # -NH2
            {'anchor': 'N', 'neighbors': ['H', 'H', 'H']}, # -NH3^+
            {'anchor': 'C', 'neighbors': ['F', 'F', 'F']},  # -CF3
            {'anchor': 'C', 'neighbors': ['H', 'H', '?']},  # -CH2?
        ]

    # Call align_with_zmat with all parameters
    align_with_zmat(
        mol1_sdf=args.molecule1,
        mol2_sdf=args.molecule2,
        processed_mol2_sdf=args.output,
        sym_equivalent_patterns=sym_equivalent_patterns,
        debug=args.debug,
        timeout=args.timeout,
        ring_matches_ring_only=args.ring_matches_ring_only,
        complete_rings_only=args.complete_rings_only,
        ring_compare=args.ring_compare,
        atom_compare=args.atom_compare,
        bond_compare=args.bond_compare,
        heavy_atom_only=args.heavy_atom_only,
        optimize_h_atoms=args.optimize_h_atoms,
        zmat_gen_method=args.zmat_gen_method,
        match_chiral_tag=args.match_chiral_tag
    )

