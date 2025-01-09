import rdkit
from rdkit import Chem
from rdkit.Chem.rdForceFieldHelpers import  MMFFOptimizeMolecule


class MDLMolecule:
    def __init__(self):
        self.atoms = []  # List of [atom_symbol, x, y, z, charge]
        self.bonds = []  # List of [atom1_idx, atom2_idx, bond_type]
        self.charges = []  # List of [atom_idx, charge]
        self.main_r_atom_index = []  # List of [atom_idx, r_group_number]  # For main structure's R# atoms
        self.attach_points = []  # List of [atom_idx, point_number]  # For R-groups
        # below is for assembled molecule
        self.main_str_atom_index = [] # List of atom indices for main structure part in assembled molecule


def find_first_ctab_section(lines, start_idx=0):
    """Find the first CTAB section starting from start_idx for main structure"""
    for i in range(start_idx, len(lines)):
        if '0999 V2000' in lines[i]:
            # Found the counts line
            counts_line = lines[i].split()
            n_atoms = int(counts_line[0])
            n_bonds = int(counts_line[1])
            return i, n_atoms, n_bonds
    return None, 0, 0


def parse_mdl_file(filename):
    """Parse MDL file into main structure and R-groups"""
    main_mol = MDLMolecule()
    r_groups = {}  # Dictionary of {r_group_num: list_of_MDLMolecule}
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    main_start, n_atoms, n_bonds = find_first_ctab_section(lines)
    if main_start is None:
        raise ValueError("No CTAB section found in file")

    # Parse main structure first
    current_line = main_start
    counts_line = lines[current_line].split()
    n_atoms = int(counts_line[0])
    n_bonds = int(counts_line[1])
    
    # Parse atoms of main structure
    current_line += 1
    for i in range(n_atoms):
        line = lines[current_line + i].split()
        atom_symbol = line[3]
        x, y, z = float(line[0]), float(line[1]), float(line[2])
        main_mol.atoms.append([atom_symbol, x, y, z, 0])
    
    # Parse bonds of main structure
    current_line += n_atoms
    for i in range(n_bonds):
        line = lines[current_line + i].split()
        atom1, atom2, bond_type = int(line[0]), int(line[1]), int(line[2])
        main_mol.bonds.append([atom1, atom2, bond_type])
    
    # Parse charges and R# atoms in main structure if present
    while current_line < len(lines) and not lines[current_line].startswith('M  END'):
        if lines[current_line].startswith('M  CHG'):
            chg_line = lines[current_line].split()
            n_charges = int(chg_line[2])
            for k in range(n_charges):
                atom_idx = int(chg_line[3+2*k])
                charge = int(chg_line[4+2*k])
                main_mol.charges.append([atom_idx, charge])

        elif lines[current_line].startswith('M  RGP'):
            # Parse R-group attachments in main structure
            line = lines[current_line]
            rgp_parts = line.split()
            n_rgroups = int(rgp_parts[2])  # Number of R-group attachments
            
            for k in range(n_rgroups):
                atom_idx = int(rgp_parts[3 + 2*k])     # Atom index in main structure
                r_group_num = int(rgp_parts[4 + 2*k])  # R-group number
                main_mol.main_r_atom_index.append([atom_idx, r_group_num])
        current_line += 1



    # Parse R-group definitions
    current_r = None
    current_r_mol = None
    
    for i, line in enumerate(lines):
        if line.startswith('$RGP'):
            # Start of a new R-group section
            current_r = int(lines[i+1].strip())
            r_groups[current_r] = []  # List to store multiple possible R-groups
            
        elif line.startswith('$CTAB') and current_r is not None:
            # Start of a new CTAB section within current R-group
            current_r_mol = MDLMolecule()
            
            # Parse atom and bond counts for this R-group
            counts_line = lines[i+1].split()  # Skip CTAB header lines
            r_n_atoms = int(counts_line[0])
            r_n_bonds = int(counts_line[1])
            
            # Parse atoms for this R-group
            for j in range(r_n_atoms):
                line = lines[i+2+j].split()
                atom_symbol = line[3]
                x, y, z = float(line[0]), float(line[1]), float(line[2])
                current_r_mol.atoms.append([atom_symbol, x, y, z, 0])
            
            # Parse bonds for this R-group
            current_line = i+2+r_n_atoms
            for j in range(r_n_bonds):
                line = lines[current_line+j].split()
                atom1, atom2, bond_type = int(line[0]), int(line[1]), int(line[2])
                current_r_mol.bonds.append([atom1, atom2, bond_type])
            
            # Parse charges and attachment points if present
            while current_line < len(lines) and not lines[current_line].startswith('M  END'):
                if lines[current_line].startswith('M  CHG'):
                    chg_line = lines[current_line].split()
                    n_charges = int(chg_line[2])
                    for k in range(n_charges):
                        atom_idx = int(chg_line[3+2*k])
                        charge = int(chg_line[4+2*k])
                        current_r_mol.charges.append([atom_idx, charge])

                elif lines[current_line].startswith('M  APO'):
                    # Parse attachment points
                    line = lines[current_line]
                    apo_parts = line.split()
                    n_attach = int(apo_parts[2])  # Number of attachment points
                    
                    for k in range(n_attach):
                        atom_idx = int(apo_parts[3 + 2*k])    # I1, I2, etc.
                        point_num = int(apo_parts[4 + 2*k])   # P1, P2, etc.
                        current_r_mol.attach_points.append([atom_idx, point_num])

                current_line += 1
            
            # Add this R-group molecule to the list for current R number
            r_groups[current_r].append(current_r_mol)
    
    return main_mol, r_groups

def assemble_molecule(main_mol, r_groups):
    """
    Assemble R groups into main structure, generating all possible combinations.
    Handles three APO scenarios:
    S1) One attach point (N=1, P1=1) - forms one bond
    S2) One attach point (N=1, P1=3) - forms two bonds
    S3) Two attach points (N=2, P1=1, P2=2) - forms two bonds with different atoms
    """
    r_group_numbers = sorted(r_groups.keys())
    variant_counts = [len(r_groups[r_num]) for r_num in r_group_numbers]
    
    import itertools
    all_combinations = []
    
    for variant_indices in itertools.product(*[range(count) for count in variant_counts]):
        # Create label for this combination
        label_parts = []
        for r_num, variant_idx in zip(r_group_numbers, variant_indices):
            label_parts.append(f"R{r_num}_{variant_idx + 1}")  # +1 for 1-based indexing
        combination_label = "_".join(label_parts)
        
        #print('molcule:', combination_label)

        variant_dict = {r_num: idx for r_num, idx in zip(r_group_numbers, variant_indices)}
        result = MDLMolecule()
        
        # Create r_assignments from main_mol.main_r_atom_index
        r_assignments = {}
        for atom_idx, r_group_num in main_mol.main_r_atom_index:
            variant_idx = variant_dict[r_group_num]
            r_assignments[atom_idx] = (r_group_num, variant_idx)
        
        atom_map = {}  # Maps old atom indices to new ones
        new_idx = 1
        
        # First, collect bonds connected to each R# atom
        r_bonds = {}  # Dictionary to store bonds for each R# atom
        for bond in main_mol.bonds:
            a1, a2, btype = bond
            if a1 in r_assignments:
                if a1 not in r_bonds:
                    r_bonds[a1] = []
                r_bonds[a1].append((a2, btype))
            elif a2 in r_assignments:
                if a2 not in r_bonds:
                    r_bonds[a2] = []
                r_bonds[a2].append((a1, btype))
        
        # Process each atom in main structure
        first_atom = True  # Flag to track first atom
        for i, atom in enumerate(main_mol.atoms, 1):
            if i in r_assignments:
                # This is an R# atom - replace with R-group
                r_num, variant_idx = r_assignments[i]
                r_group = r_groups[r_num][variant_idx]
                
                # Get attachment points
                attach_points = r_group.attach_points
                if not attach_points:
                    raise ValueError(f"No attachment points found in R-group {r_num}")
                
                connected_bonds = r_bonds.get(i, [])
                
                # Handle different APO scenarios
                if len(attach_points) == 1:
                    # S1 or S2
                    attach_idx, point_num = attach_points[0]
                    if point_num == 1:
                        # S1: One attachment point, one bond
                        if len(connected_bonds) != 1:
                            raise ValueError(f"Scenario S1: Expected 1 bond, found {len(connected_bonds)} for R-group {r_num}")
                    elif point_num == 3:
                        # S2: One attachment point, two bonds
                        if len(connected_bonds) != 2:
                            raise ValueError(f"Scenario S2: Expected 2 bonds, found {len(connected_bonds)} for R-group {r_num}")
                elif len(attach_points) == 2:
                    # S3: Two attachment points, two bonds
                    if len(connected_bonds) != 2:
                        raise ValueError(f"Scenario S3: Expected 2 bonds, found {len(connected_bonds)} for R-group {r_num}")
                    # Sort attachment points by point number (1 and 2)
                    attach_points = sorted(attach_points, key=lambda x: x[1])
                
                # Add all atoms from R-group
                for j, r_atom in enumerate(r_group.atoms, 1):
                    if first_atom:
                        # Add 0.01 to z-coordinate of first atom
                        modified_atom = r_atom.copy()  # Create a copy of the atom list
                        modified_atom[3] += 0.01  # Increment z-coordinate
                        result.atoms.append(modified_atom)
                        first_atom = False
                    else:
                        result.atoms.append(r_atom)
                    atom_map[f"R{r_num}_{j}"] = new_idx
                    new_idx += 1
                
                # Create bonds between attachment points and main structure
                if len(attach_points) == 1 and attach_points[0][1] == 3:
                    # S2: Both bonds connect to the same attachment point
                    attach_idx = attach_points[0][0]
                    new_attach_idx = atom_map[f"R{r_num}_{attach_idx}"]
                    for other_atom, bond_type in connected_bonds:
                        new_other_idx = atom_map.get(other_atom, None)
                        if new_other_idx is None:
                            atom_map[other_atom] = new_idx
                            result.atoms.append(main_mol.atoms[other_atom-1])
                            new_other_idx = new_idx
                            result.main_str_atom_index.append(new_idx)
                            new_idx += 1
                        result.bonds.append([new_attach_idx, new_other_idx, bond_type])
                else:
                    # S1 or S3: One bond per attachment point
                    for (attach_idx, point_num), (other_atom, bond_type) in zip(attach_points, connected_bonds):
                        new_attach_idx = atom_map[f"R{r_num}_{attach_idx}"]
                        new_other_idx = atom_map.get(other_atom, None)
                        if new_other_idx is None:
                            atom_map[other_atom] = new_idx
                            result.atoms.append(main_mol.atoms[other_atom-1])
                            new_other_idx = new_idx
                            result.main_str_atom_index.append(new_idx)
                            new_idx += 1
                        result.bonds.append([new_attach_idx, new_other_idx, bond_type])
                
                # Add internal R-group bonds
                for bond in r_group.bonds:
                    a1, a2, btype = bond
                    new_a1 = atom_map[f"R{r_num}_{a1}"]
                    new_a2 = atom_map[f"R{r_num}_{a2}"]
                    result.bonds.append([new_a1, new_a2, btype])
                
                # Add R-group charges
                for charge_idx, charge in r_group.charges:
                    new_charge_idx = atom_map[f"R{r_num}_{charge_idx}"]
                    result.charges.append([new_charge_idx, charge])
                    
            else:
                # Regular atom - copy as is
                if i not in atom_map:
                    if first_atom:
                        # Add 0.01 to z-coordinate of first atom
                        modified_atom = atom.copy()  # Create a copy of the atom list
                        modified_atom[3] += 0.01  # Increment z-coordinate
                        result.atoms.append(modified_atom)
                        first_atom = False
                    else:
                        result.atoms.append(atom)
                    atom_map[i] = new_idx
                    result.main_str_atom_index.append(new_idx)
                    new_idx += 1
        
        # Copy remaining main structure bonds
        for bond in main_mol.bonds:
            a1, a2, btype = bond
            if a1 not in r_assignments and a2 not in r_assignments:
                new_a1 = atom_map[a1]
                new_a2 = atom_map[a2]
                result.bonds.append([new_a1, new_a2, btype])
        
        # Copy main structure charges
        for charge_idx, charge in main_mol.charges:
            if charge_idx not in r_assignments:
                new_charge_idx = atom_map[charge_idx]
                result.charges.append([new_charge_idx, charge])
        
        all_combinations.append((result, combination_label))
    
    return all_combinations


def rdkit_optimize_molecule(input_sdf, output_sdf, optimize_h_only=False):
    """
    Optimize molecule using RDKit with MMFF force field and fix charge indicators
    Args:
        input_sdf: input SDF file
        output_sdf: output SDF file    
        optimize_h_only: if True, only optimize H atoms 
    """
    # First let RDKit optimize and write the molecule
    mm = rdkit.Chem.SDMolSupplier(input_sdf, removeHs=False, sanitize=True)
    mol = mm[0]

    if optimize_h_only:
        mp = Chem.AllChem.MMFFGetMoleculeProperties(mol)
        ff = Chem.AllChem.MMFFGetMoleculeForceField(mol, mp)
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != 1:  # If not hydrogen
                ff.AddFixedPoint(atom.GetIdx())
        ff.Minimize(maxIts=200)
    else:
        MMFFOptimizeMolecule(mol)


    w = Chem.SDWriter(output_sdf)
    w.write(mol)
    w.close()

    # Now read and modify the generated SDF file
    charge_to_indicator = {
        -3: 7, -2: 6, -1: 5, 0: 0, 1: 3, 2: 2, 3: 1
    }
    
    # Read the file and store all lines
    with open(output_sdf, 'r') as f:
        lines = f.readlines()

    # Find charges from 'M  CHG' lines
    atom_charges = {}
    for line in lines:
        if line.startswith('M  CHG'):
            parts = line.split()
            n_charges = int(parts[2])
            for i in range(n_charges):
                atom_idx = int(parts[3 + 2*i])
                charge = int(parts[4 + 2*i])
                atom_charges[atom_idx] = charge

    # Find the atom block and modify charge indicators
    for i, line in enumerate(lines):
        if 'V2000' in line:
            atom_block_start = i + 1
            n_atoms = int(line.split()[0])
            # Modify each atom line
            for j in range(n_atoms):
                line_idx = atom_block_start + j
                atom_line = lines[line_idx]
                # Split the line maintaining the exact formatting
                coords_and_symbol = atom_line[:36]  # coordinates and atom symbol plus first number
                charge_indicator = charge_to_indicator.get(atom_charges.get(j + 1, 0), 0)
                rest_of_line = atom_line[39:]  # rest of the line after charge indicator
                # Reconstruct the line with the new charge indicator
                lines[line_idx] = f"{coords_and_symbol}{charge_indicator:3d}{rest_of_line}"
            break

    # Write the modified file
    with open(output_sdf, 'w') as f:
        f.writelines(lines)

    return


# file I/O

def write_mol_file(molecule, filename, label=None):
    """Write assembled molecule to MOL file"""
    # Create charge indicator mapping
    charge_to_indicator = {
        -3: 7,
        -2: 6,
        -1: 5,
        0: 0,
        1: 3,
        2: 2,
        3: 1
    }
    
    # Create a lookup dictionary for atom charges
    atom_charges = {atom_idx: charge for atom_idx, charge in molecule.charges}
    
    with open(filename, 'w') as f:
        # Write header
        if label:
            f.write(label+"\n\n\n")
        else:
            f.write("assembled molecule\n\n\n")
        
        # Write counts line
        n_atoms = len(molecule.atoms)
        n_bonds = len(molecule.bonds)
        f.write(f"{n_atoms:3d}{n_bonds:3d}  0  0  0  0  0  0  0  0999 V2000\n")
        
        # Write atoms
        for i, atom in enumerate(molecule.atoms, 1):
            symbol, x, y, z, _ = atom
            # Get charge indicator (default to 0 if no charge)
            charge = atom_charges.get(i, 0)
            charge_indicator = charge_to_indicator.get(charge, 0)
            f.write(f"{x:10.4f}{y:10.4f}{z:10.4f} {symbol:<3s} 0  {charge_indicator}  0  0  0  0  0  0  0  0  0  0\n")
        
        # Write bonds
        for bond in molecule.bonds:
            a1, a2, btype = bond
            f.write(f"{a1:3d}{a2:3d}{btype:3d}  0  0  0  0\n")
        
        # Write charges if present
        if molecule.charges:
            # Format: M  CHG  n atom1 charge1 atom2 charge2 ... atomn chargen
            # where n is the number of charged atoms
            n_charges = len(molecule.charges)
            charge_str = f"M  CHG  {n_charges}"
            for atom_idx, charge in molecule.charges:
                charge_str += f" {atom_idx:3d} {charge:3d}"
            f.write(charge_str + "\n")
        
        # Write M END
        f.write("M  END\n")
 

