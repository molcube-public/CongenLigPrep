import os
import argparse
import sys
import tempfile
import shutil
import json

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(script_dir)
from core_3 import process_smiles_file
from core_1 import align_with_zmat, str2bool


def process_and_align_smiles(input_smi, align=False, reference_sdf=None, output_base='molecule',
                           timeout=10, ring_matches_ring_only=True,
                           complete_rings_only=False, ring_compare='IgnoreRingFusion',
                           atom_compare='CompareElements', bond_compare='CompareOrder',
                           debug=0, intermedia_output=None, sym_equivalent_patterns=None,
                           heavy_atom_only=False, optimize_h_atoms=False, zmat_gen_method=1, match_chiral_tag=True):
    """
    Process SMILES file to generate 3D SDF files and optionally align them to a reference structure.
    
    Args:
        input_smi (str): Path to input SMILES file
        align (bool): Whether to perform alignment
        reference_sdf (str, optional): Path to reference SDF file for alignment
        output_base (str): Base name for output files
        timeout (int): Maximum time in seconds for MCS search
        ring_matches_ring_only (bool): If True, ring atoms only match ring atoms
        complete_rings_only (bool): If True, only complete rings are considered
        ring_compare (str): Ring matching behavior
        atom_compare (str): Atom comparison method
        bond_compare (str): Bond comparison method
        debug (int): Debug level (0-2)
        intermedia_output (str, optional): Path to save intermediate SDF files
        sym_equivalent_patterns (list, optional): List of patterns for symmetry equivalent groups
            e.g. [{'anchor': 'C', 'neighbors': ['H', 'H', 'H']}, ...]
        heavy_atom_only (bool): If True, only heavy atoms are considered during MCS search
        optimize_h_atoms (bool): If True, optimize H atoms when heavy_atom_only is True
        match_chiral_tag (bool): If True, chiral tags are matched
    Returns:
        list: Paths to generated SDF files
    """
    # Initialize temp_dir
    temp_dir = None
    
    if align:
        # Only create temp directory or use intermedia_output when aligning
        if intermedia_output is None:
            temp_dir = tempfile.mkdtemp()
            working_output_base = os.path.join(temp_dir, os.path.basename(output_base))
        else:
            # Use intermedia_output directly as the base name
            os.makedirs(os.path.dirname(intermedia_output) or '.', exist_ok=True)
            working_output_base = intermedia_output
    else:
        # If not aligning, output directly to final destination
        working_output_base = output_base
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_base) or '.', exist_ok=True)
    
    try:
        # Generate SDF files from SMILES
        sdf_files = process_smiles_file(input_smi, working_output_base)
        print(len(sdf_files), 'SMILES entries have been converted to SDF files.')
        
        # Convert relative paths to absolute paths
        sdf_files = [os.path.abspath(f) for f in sdf_files]
        
        if not sdf_files:
            raise RuntimeError("No SDF files were generated from the SMILES file")
            
        if not align:
            return sdf_files  # Return directly as files are already in correct location
            
        # Determine reference structure
        if reference_sdf is None:
            reference_sdf = sdf_files[0]
            files_to_align = sdf_files[1:]
            print(f"\nUsing {reference_sdf} as reference structure")
        else:
            files_to_align = sdf_files  # All files need to be aligned when using external reference
            print(f"\nUsing external reference structure: {reference_sdf}")
        
        aligned_files = []
        for i, sdf_file in enumerate(files_to_align):
            # Get the molecule number, accounting for reference structure
            if reference_sdf == sdf_files[0]:
                mol_num = i + 2  # Start from 2 when using first molecule as reference
            else:
                mol_num = i + 1  # Start from 1 when using external reference
            
            # Create output filename using output_base path
            aligned_name = f"{os.path.basename(output_base)}_{mol_num}_aligned.sdf"
            aligned_sdf = os.path.join(os.path.dirname(output_base), aligned_name)
            
            # Align the structure
            try:
                align_with_zmat(
                    mol1_sdf=reference_sdf,
                    mol2_sdf=sdf_file,
                    processed_mol2_sdf=aligned_sdf,
                    sym_equivalent_patterns=sym_equivalent_patterns,
                    debug=debug,
                    timeout=timeout,
                    ring_matches_ring_only=ring_matches_ring_only,
                    complete_rings_only=complete_rings_only,
                    ring_compare=ring_compare,
                    atom_compare=atom_compare,
                    bond_compare=bond_compare,
                    heavy_atom_only=heavy_atom_only,
                    optimize_h_atoms=optimize_h_atoms,
                    zmat_gen_method=zmat_gen_method,
                    match_chiral_tag=match_chiral_tag
                )
                aligned_files.append(aligned_sdf)
            except Exception as e:
                print(f"Warning: Failed to align {sdf_file}: {str(e)}")
                # Copy original file to output directory if alignment fails
                final_path = os.path.join(os.path.dirname(output_base), 
                                        f"{os.path.basename(output_base)}_{mol_num}.sdf")
                shutil.copy2(sdf_file, final_path)
                aligned_files.append(final_path)
        
        # If using internal reference (first structure), copy it to output directory
        if reference_sdf == sdf_files[0]:
            final_ref = os.path.join(os.path.dirname(output_base), 
                                   f"{os.path.basename(output_base)}_1.sdf")
            shutil.copy2(reference_sdf, final_ref)
            return [final_ref] + aligned_files
        
        return aligned_files
        
    finally:
        # Clean up temp directory if created
        if temp_dir:
            shutil.rmtree(temp_dir)

if __name__ == '__main__':
    #parser = argparse.ArgumentParser(description='Convert SMILES to 3D SDF files with optional alignment')
    class ArgumentParserWithValidation(argparse.ArgumentParser):
        def parse_args(self):
            args = super().parse_args()
            if args.optimize_h_atoms and not args.heavy_atom_only:
                self.error("--optimize-h-atoms can only be used when --heavy-atom-only is True")
            return args

    parser = ArgumentParserWithValidation(description='Convert SMILES to 3D SDF files with optional alignment')
    
    # Required arguments
    parser.add_argument('-i', '--input', required=True,
                      help='Input SMILES file path')
    
    # Optional arguments
    parser.add_argument('-o', '--output',
                      help='Output base name for SDF files (defaults to input filename without extension)')
    parser.add_argument('--align', action='store_true',
                      help='Align structures (if not specified, only generate 3D structures)')
    parser.add_argument('-r', '--reference',
                      help='Reference SDF file for alignment (if not specified but --align is used, first generated structure will be reference)')
    
    # Alignment options
    align_group = parser.add_argument_group('Alignment options (only used with --align)')
    align_group.add_argument('--timeout', type=int, default=10,
                      help='Maximum time in seconds for MCS search (default: 10)')
    align_group.add_argument('--ring-matches-ring-only', type=str2bool, default=True, metavar='BOOL',
                      help='Ring atoms only match ring atoms (default: True)')
    align_group.add_argument('--complete-rings-only', type=str2bool, default=False, metavar='BOOL',
                      help='Only complete rings are considered (default: False)')
    align_group.add_argument('--ring-compare', 
                      choices=['StrictRingFusion', 'IgnoreRingFusion', 'PermissiveRingFusion'],
                      default='IgnoreRingFusion',
                      help='Ring matching behavior (default: IgnoreRingFusion)')
    align_group.add_argument('--atom-compare',
                      choices=['CompareAny', 'CompareElements', 'CompareIsotopes', 'CompareAnyHeavyAtom'],
                      default='CompareElements',
                      help='Atom comparison method (default: CompareElements)')
    align_group.add_argument('--bond-compare',
                      choices=['CompareAny', 'CompareOrder', 'CompareOrderExact'],
                      default='CompareOrder',
                      help='Bond comparison method (default: CompareOrder)')
    align_group.add_argument('--debug', '-d', type=int, default=0,
                      help='Debug level (0: disabled, 1: basic info, 2: basic + 2D molecule visualization) (default: 0)')
    
    # Add new argument for intermediate files, only used with --align
    align_group.add_argument('--intermedia-output',
                      help='Base name for intermediate SDF files (only used with --align, if not specified, uses temp directory)')
    
    # Add patterns argument to alignment options group
    align_group.add_argument('--patterns', '-p',
                      help='Path to JSON file containing symmetry equivalent patterns')
    
    # Add heavy-atom-only argument to alignment options group
    align_group.add_argument('--heavy-atom-only', type=str2bool, default=False, metavar='BOOL',
                      help='Only heavy atoms are considered during MCS search (default: False)')
    
    # Add optimize-h-atoms argument to alignment options group
    align_group.add_argument('--optimize-h-atoms', type=str2bool, default=False, metavar='BOOL',
                      help='Optimize H atoms (default: False)')

    # Add zmat-gen-method argument to alignment options group
    align_group.add_argument('--zmat-gen-method', 
                      choices=[1, 2, 3], 
                      type=int,
                      default=1,
                      help='Method to generate z-matrix (1: use our own zmat generation, 2: use geomConvert\'s method, 3: use XYZ-to-ZMAT\'s method) (default: 1)')
    
    # Add match-chiral-tag argument to alignment options group
    align_group.add_argument('--match-chiral-tag', type=str2bool, default=True, metavar='BOOL',
                      help='Match chiral tags (default: True)')
    
    args = parser.parse_args()
    
    # Get default output name from input filename if not specified
    output_base = args.output if args.output else os.path.splitext(os.path.basename(args.input))[0]
    
    # Load patterns from file if provided, otherwise use defaults
    sym_equivalent_patterns = None
    if args.align:  # Only process patterns if alignment is requested
        if args.patterns:
            with open(args.patterns, 'r') as f:
                sym_equivalent_patterns = json.load(f)
    
    # Process SMILES and align if requested
    output_files = process_and_align_smiles(
        input_smi=args.input,
        align=args.align,
        reference_sdf=args.reference,
        output_base=output_base,
        timeout=args.timeout,
        ring_matches_ring_only=args.ring_matches_ring_only,
        complete_rings_only=args.complete_rings_only,
        ring_compare=args.ring_compare,
        atom_compare=args.atom_compare,
        bond_compare=args.bond_compare,
        debug=args.debug,
        intermedia_output=args.intermedia_output,
        sym_equivalent_patterns=sym_equivalent_patterns,
        heavy_atom_only=args.heavy_atom_only,
        optimize_h_atoms=args.optimize_h_atoms,
        zmat_gen_method=args.zmat_gen_method,
        match_chiral_tag=args.match_chiral_tag
    )
    
    # Print output files
    print("\nGenerated SDF files:")
    for f in output_files:
        print(f)

