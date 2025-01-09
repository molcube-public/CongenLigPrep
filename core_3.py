import os
import re
import argparse

def process_smiles_file(input_smi, output_base='molecule'):
    # Read the original SMILES file and store molecule labels if they exist
    # convert SMILES to SDF files using OpenBabel
    
    mol_labels = []
    cleaned_smiles = []
    
    with open(input_smi, 'r') as f:
        for line in f:
            parts = line.strip().split(',')
            if len(parts) >= 2:  # If there's a label (and possibly more fields)
                smiles_part, label = parts[0:2]  # Take only first two fields
                mol_labels.append(label)
                cleaned_smiles.append(smiles_part)
            else:  # If there's only SMILES
                cleaned_smiles.append(parts[0])
    
    sdf_files = []
    # Process each SMILES string individually, using idx+1 for filenames
    for idx, smiles in enumerate(cleaned_smiles):
        # Create individual temporary SMILES file with 1-based indexing
        temp_smi = f'temp_cleaned_{idx+1}.smi'
        output_sdf = f'{output_base}_{idx+1}.sdf'
        
        try:
            # Write single SMILES to temporary file
            with open(temp_smi, 'w') as f:
                f.write(f"{smiles}\n")
            
            # Run OpenBabel for single molecule
            cmd = f'obabel -i smi {temp_smi} -o sdf -O {output_sdf} --gen3d slow > /dev/null 2>&1'
            if os.system(cmd) != 0:
                raise RuntimeError(f"OpenBabel command failed for SMILES entry no. {idx+1}")
            
            # Process labels if they exist
            if mol_labels:
                with open(output_sdf, 'r') as f:
                    content = f.read()
                modified_content = re.sub(r'^.*?\n', f"{mol_labels[idx]}\n", content, count=1)
                with open(output_sdf, 'w') as f:
                    f.write(modified_content)
            
            sdf_files.append(output_sdf)
            
        finally:
            # Clean up temporary file
            if os.path.exists(temp_smi):
                os.remove(temp_smi)

    return sdf_files

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process SMILES files to generate 3D SDF files.')
    parser.add_argument('-i', '--input', required=True,
                      help='Input SMILES file path')
    parser.add_argument('-o', '--output', 
                      help='Output base name for SDF files (defaults to input filename without extension)')
    
    args = parser.parse_args()
    # Get default output name from input filename if not specified
    output_base = args.output if args.output else os.path.splitext(os.path.basename(args.input))[0]
    process_smiles_file(args.input, output_base) 

