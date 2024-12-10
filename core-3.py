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
            if len(parts) > 1:  # If there's a label
                smiles_part, label = parts
                mol_labels.append(label)
                cleaned_smiles.append(smiles_part)
            else:  # If there's only SMILES
                cleaned_smiles.append(parts[0])
    
    # Write temporary SMILES file without labels
    temp_smi = 'temp_cleaned.smi'
    with open(temp_smi, 'w') as f:
        for smiles in cleaned_smiles:
            f.write(f"{smiles}\n")
    
    # Check if obabel exists and works
    if os.system('which obabel >/dev/null 2>&1') != 0:
        raise RuntimeError("OpenBabel (obabel) is not found in PATH. Please install OpenBabel or check your PATH settings.")

    # Run OpenBabel
    cmd = f'obabel -i smi {temp_smi} -o sdf -O {output_base}_.sdf --gen3d -m'
    if os.system(cmd) != 0:
        os.remove(temp_smi)  # Clean up temp file if OpenBabel fails
        raise RuntimeError("OpenBabel command failed")
    
    # Get list of generated SDF files
    def get_number(filename):
        match = re.search(r'.*_(\d+)\.sdf$', filename)
        return int(match.group(1)) if match else 0
        
    sdf_files = [f for f in os.listdir('.') if f.startswith(f'{output_base}_') and f.endswith('.sdf')]
    sdf_files.sort(key=get_number)
    
    # Process each generated SDF file only if we have labels
    if mol_labels:
        for i, sdf_file in enumerate(sdf_files):
            with open(sdf_file, 'r') as f:
                content = f.read()
            
            # Replace the first line with the molecule label
            modified_content = re.sub(r'^.*?\n', f"{mol_labels[i]}\n", content, count=1)
            
            with open(sdf_file, 'w') as f:
                f.write(modified_content)
    
    # Clean up temporary file
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

