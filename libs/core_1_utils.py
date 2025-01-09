import numpy as np
from scipy.sparse.csgraph import connected_components
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS
import os
import tempfile

from gcutil import *
from xyz_to_zmat import xyz_to_zmat 

def find_mcs_matches_rdkit_2d(mol1, mol2, timeout=10, 
                     ring_matches_ring_only=True,
                     complete_rings_only=False,
                     ring_compare='IgnoreRingFusion',
                     atom_compare='CompareElements',
                     bond_compare='CompareOrderExact',
                     match_chiral_tag=True):
    """
    Find Maximum Common Substructure (MCS) between two molecules and return matching atom indices.
    
    Args:
        mol1 (rdkit.Chem.rdchem.Mol): First molecule
        mol2 (rdkit.Chem.rdchem.Mol): Second molecule
        timeout (int): Maximum time in seconds for MCS search
        ring_matches_ring_only (bool): If True, ring atoms are only allowed to match ring atoms
        complete_rings_only (bool): If True, only complete rings are considered in matches
        ring_compare (str): Ring matching behavior. Options:
            - 'StrictRingFusion'
            - 'IgnoreRingFusion'
            - 'PermissiveRingFusion'
        atom_compare (str): Atom comparison method. Options:
            - 'CompareAny'
            - 'CompareElements'
            - 'CompareIsotopes'
            - 'CompareAnyHeavyAtom'
        bond_compare (str): Bond comparison method. Options:
            - 'CompareAny'
            - 'CompareOrder'
            - 'CompareOrderExact'
        match_chiral_tag (bool): If True, chiral tags are matched
    Returns:
        list: List of tuples containing matching atom indices for each molecule
    """
    # Define mappings from strings to RDKit enums
    ring_compare_options = {
        'StrictRingFusion': rdkit.Chem.rdFMCS.RingCompare.StrictRingFusion,
        'IgnoreRingFusion': rdkit.Chem.rdFMCS.RingCompare.IgnoreRingFusion,
        'PermissiveRingFusion': rdkit.Chem.rdFMCS.RingCompare.PermissiveRingFusion
    }
    
    atom_compare_options = {
        'CompareAny': rdkit.Chem.rdFMCS.AtomCompare.CompareAny,
        'CompareAnyHeavyAtom':rdkit.Chem.rdFMCS.AtomCompare.CompareAnyHeavyAtom,
        'CompareElements': rdkit.Chem.rdFMCS.AtomCompare.CompareElements,
        'CompareIsotopes': rdkit.Chem.rdFMCS.AtomCompare.CompareIsotopes
    }
    
    bond_compare_options = {
        'CompareAny': rdkit.Chem.rdFMCS.BondCompare.CompareAny,
        'CompareOrder': rdkit.Chem.rdFMCS.BondCompare.CompareOrder,
        'CompareOrderExact': rdkit.Chem.rdFMCS.BondCompare.CompareOrderExact
    }
    
    # Convert string inputs to RDKit enums
    try:
        ring_compare_enum = ring_compare_options[ring_compare]
        atom_compare_enum = atom_compare_options[atom_compare]
        bond_compare_enum = bond_compare_options[bond_compare]
    except KeyError as e:
        raise ValueError(f"Invalid option provided: {e}. Please check the documentation for valid options.")
    
    fmcs_params = dict(
        timeout=timeout,
        ringMatchesRingOnly=ring_matches_ring_only,
        completeRingsOnly=complete_rings_only,
        ringCompare=ring_compare_enum,
        atomCompare=atom_compare_enum,
        bondCompare=bond_compare_enum,
        matchValences=True, # newly added 
        matchChiralTag=match_chiral_tag, # newly added 
    )
    
    r = rdkit.Chem.rdFMCS.FindMCS([mol1, mol2], **fmcs_params)
    patt = Chem.MolFromSmarts(r.smartsString)
    
    hits = []
    hits.append(mol1.GetSubstructMatch(patt))
    hits.append(mol2.GetSubstructMatch(patt))
    
    return hits, r


def find_substituent_indices(mol, mcs_hits):
    """
    Find substituent indices in a molecule by analyzing its connectivity relative to MCS.
    
    Args:
        mol (rdkit.Chem.rdchem.Mol): Input RDKit molecule
        mcs_hits (list): List of atom indices that are part of MCS
    
    Returns:
        dict: Dictionary where keys are substituent site labels and values are lists of atom indices
    """
    # Get adjacency matrix
    adj_matrix = Chem.GetAdjacencyMatrix(mol)
    
    # Erase MCS connectivity info
    setZeroRemoveCoreExtend(adj_matrix, mcs_hits)
    
    # Find connected components
    num_components, labels_scatter = connected_components(adj_matrix)
    
    # Remove MCS atoms from labels
    labels_scatter = rmMCSfromConnComp(mcs_hits, labels_scatter)
    
    # Get unique labels excluding MCS (-1)
    uniq_lab = list(set(labels_scatter))
    uniq_lab.remove(-1)
    
    # Create dictionary of substituent indices
    substituent_indices = {}
    for label in uniq_lab:
        indices = np.where(labels_scatter == label)[0]
        substituent_indices[label] = indices.tolist()
    
    return substituent_indices


def find_anchor_and_neighbors(mol, substituent_atoms, mcs_atoms, n_total=3):
    """
    For each substituent, find anchor atoms (MCS atoms directly connected to substituent)
    and their closest neighbors in the MCS to get a total of n_total atoms.
    
    Args:
        mol: RDKit molecule
        substituent_atoms: list of atom indices for the substituent
        mcs_atoms: list of atom indices for the MCS
        n_total: total number of atoms to return (anchor + neighbors)
    
    Returns:
        list of atom indices [anchor_atoms + neighbor_atoms]
    """
    # Get adjacency matrix
    adj_matrix = Chem.GetAdjacencyMatrix(mol)
    
    # Find first anchor atom (MCS atom directly connected to substituent)
    anchor_atoms = []
    found_anchor = False
    for mcs_idx in mcs_atoms:
        if found_anchor:
            break
        for sub_idx in substituent_atoms:
            if adj_matrix[mcs_idx][sub_idx] == 1:
                anchor_atoms.append(mcs_idx)
                found_anchor = True
                break
    
    # Get 3D coordinates for distance calculation
    conf = mol.GetConformer()
    pos = conf.GetPositions()
    
    # Find neighboring atoms in MCS
    n_neighbors_needed = n_total - len(anchor_atoms)
    if n_neighbors_needed > 0:
        # Calculate distances from anchor atoms to all MCS atoms
        distances = []
        potential_neighbors = [x for x in mcs_atoms if x not in anchor_atoms]
        
        for mcs_idx in potential_neighbors:
            min_dist = float('inf')
            for anchor_idx in anchor_atoms:
                dist = np.linalg.norm(pos[mcs_idx] - pos[anchor_idx])
                min_dist = min(min_dist, dist)
            distances.append((mcs_idx, min_dist))
        
        # Sort by distance and get the closest neighbors
        distances.sort(key=lambda x: x[1])
        neighbor_atoms = [x[0] for x in distances[:n_neighbors_needed]]
    else:
        neighbor_atoms = []
    
    return anchor_atoms + neighbor_atoms



def collect_substituent_info(mol, mol_2d, substituent_indices, mcs_atoms, visualize=False, debug=False):
    """
    Analyze substituents by finding anchor atoms and their neighbors for each substituent site.
    
    Args:
        mol: RDKit 3D molecule
        mol_2d: RDKit 2D molecule for visualization
        substituent_indices: dictionary of substituent sites and their atom indices
        mcs_atoms: list of MCS atom indices
        visualize: boolean to control visualization output
        debug: boolean to control debug print statements
    
    Returns:
        dictionary containing substituent analysis results
    """
    results = {}
    for site_label, sub_atoms in substituent_indices.items():
        anchor_and_neighbors = find_anchor_and_neighbors(mol, sub_atoms, mcs_atoms)
        results[site_label] = {
            'substituent_atoms': sub_atoms,
            'anchor_and_neighbors': anchor_and_neighbors
        }
        if debug:
            print(f"\nSubstituent site {site_label}:")
            print(f"  Substituent atoms: {sub_atoms}")
            print(f"  Anchor + neighbor atoms: {anchor_and_neighbors}")

    if visualize:
        for site_label, data in results.items():
            display('sub+anc+neigh:')
            display(rdkit.Chem.Draw.MolToImage(
                mol_2d, 
                highlightAtoms=data['substituent_atoms'] + data['anchor_and_neighbors']
            ))
            display('sub:')
            display(rdkit.Chem.Draw.MolToImage(
                mol_2d, 
                highlightAtoms=data['substituent_atoms'] 
            ))
    
    return results


def order_atoms_by_distance(elements, coordinates, property=None):
    """
    Reorder atoms based on spatial proximity, keeping the first atom fixed.
    
    Args:
        elements (list): List of element symbols
        coordinates (list): List of [x,y,z] coordinates for each atom
        property (list, optional): List of properties for each atom
    
    Returns:
        tuple: (ordered_elements, ordered_coordinates, ordered_properties)
            - ordered_elements: List of reordered element symbols
            - ordered_coordinates: List of reordered coordinates
            - ordered_properties: List of reordered properties (if provided)
    """
    
    
    coords = np.array(coordinates)
    n_atoms = len(elements)
    
    # Initialize ordered indices with the first atom
    ordered_indices = [0]
    remaining_indices = set(range(1, n_atoms))
    
    # Process atoms until all are ordered
    while remaining_indices:
        ref_idx = ordered_indices[-1]
        ref_pos = coords[ref_idx]
        
        # Calculate distances from reference atom to all remaining atoms
        distances = []
        for idx in remaining_indices:
            dist = np.linalg.norm(coords[idx] - ref_pos)
            distances.append((idx, dist))
        
        # Sort distances
        distances.sort(key=lambda x: x[1])
        
        # First handle atoms within 1 Å (if any)
        close_atoms = [idx for idx, dist in distances if dist < 1.0]
        for idx in close_atoms:
            ordered_indices.append(idx)
            remaining_indices.remove(idx)
            
        # If no close atoms found (or after processing them),
        # add the next closest atom
        if remaining_indices:
            for idx, dist in distances:
                if idx in remaining_indices:
                    ordered_indices.append(idx)
                    remaining_indices.remove(idx)
                    break
    
    # Reorder all data using the new indices
    ordered_elements = [elements[i] for i in ordered_indices]
    ordered_coordinates = [coordinates[i] for i in ordered_indices]
    
    if property is not None:
        ordered_properties = [property[i] for i in ordered_indices]
        return ordered_elements, ordered_coordinates, ordered_properties
    
    return ordered_elements, ordered_coordinates, None




def extract_anchor_subst_xyz(mol, substituent_data):
    """
    Collect atom elements and XYZ coordinates for a substituent site and its anchor/neighbor atoms.
    
    Args:
        mol: RDKit 3D molecule
        substituent_data: dictionary containing 'substituent_atoms' and 'anchor_and_neighbors'
            from collect_substituent_info results
    
    Returns:
        dict: Dictionary containing:
            - 'elements': list of atom elements in order [anchor, neighbors, substituent]
            - 'coordinates': numpy array of XYZ coordinates in same order
            - 'indices': list of atom indices [anchor, neighbors, substituent]
    """

    debug = 0

    # Split anchor_and_neighbors into anchor atoms and neighbor atoms
    # (anchor atoms are those directly connected to substituent)
    adj_matrix = Chem.GetAdjacencyMatrix(mol)
    anchor_and_neighbors = substituent_data['anchor_and_neighbors']
    substituent_atoms = substituent_data['substituent_atoms']
    
    # Identify anchor atoms (those connected to substituent)
    assert len(anchor_and_neighbors) == 3
    
    anchor_atoms = []
    neighbor_atoms = []
    for atom_idx in anchor_and_neighbors:
        # Skip if we already found an anchor atom
        if len(anchor_atoms) >= 1:
            neighbor_atoms.append(atom_idx)
            continue
            
        is_anchor = False
        for sub_idx in substituent_atoms:
            if adj_matrix[atom_idx][sub_idx] == 1:
                is_anchor = True
                break
        if is_anchor:
            anchor_atoms.append(atom_idx)
        else:
            neighbor_atoms.append(atom_idx)
    

    assert len(anchor_atoms) == 1
    # Combine all indices in desired order
    ordered_indices = anchor_atoms + neighbor_atoms + substituent_atoms
    
    # Get conformer for 3D coordinates
    conf = mol.GetConformer()
    
    # Collect elements and coordinates
    elements_wo_ordering = []
    coordinates_wo_ordering = []
    


    for idx in ordered_indices:
        atom = mol.GetAtomWithIdx(idx)
        elements_wo_ordering.append(atom.GetSymbol())
        pos = conf.GetAtomPosition(idx)
        coordinates_wo_ordering.append([pos.x, pos.y, pos.z])

    # debug print:
    if debug:
        print('site_xyz before ordering:')
        print(len(elements_wo_ordering))
        print('')
        for i in range(len(elements_wo_ordering)):
            print(elements_wo_ordering[i], coordinates_wo_ordering[i][0], coordinates_wo_ordering[i][1], coordinates_wo_ordering[i][2])
        print('end.')


    # prepare for ordering by distance
    indices_to_send = anchor_atoms + substituent_atoms
    elements_to_send = []
    coordinates_to_send = []
    for idx in indices_to_send:
        atom = mol.GetAtomWithIdx(idx)
        elements_to_send.append(atom.GetSymbol())
        pos = conf.GetAtomPosition(idx)
        coordinates_to_send.append([pos.x, pos.y, pos.z])

    # order by distance
    elements_anc_sub, coordinates_anc_sub, indices_anc_sub = order_atoms_by_distance(elements_to_send, coordinates_to_send, indices_to_send)

    # combine things together 
    elements = elements_wo_ordering[:3] + elements_anc_sub[1:]
    coordinates = coordinates_wo_ordering[:3] + coordinates_anc_sub[1:]  
    indices = ordered_indices[:3] + indices_anc_sub[1:]


    
    # debug print:
    if debug:
        print('site_xyz:')
        print(len(elements))
        print('')
        for i in range(len(elements)):
            print(elements[i], coordinates[i][0], coordinates[i][1], coordinates[i][2])
        print('end.')
    
    


    return {
        'elements': elements,
        'coordinates': np.array(coordinates),
        'indices': indices
    }




def update_coordinates(input_coords, indices, new_coordinates):
    """
    Update specific coordinates in the original coordinate array.
    
    Args:
        orig_coords (numpy.ndarray): Original coordinates array of shape (n_atoms, 3)
        indices (list): List of atom indices to modify
        new_coordinates (numpy.ndarray): New coordinates array of shape (n_indices, 3)
        
    Returns:
        numpy.ndarray: Updated coordinates array with same shape as orig_coords
    """
    ## Make a copy of original coordinates to avoid modifying the input
    #new_coords = orig_coords.copy()
    
    # Update coordinates at specified indices
    for i, idx in enumerate(indices):
        input_coords[idx] = new_coordinates[i]
    
    return input_coords



def strip_sym_equiv_hydrogens(mol, hits, mol_index):
    """
    Remove symmetrically equivalent hydrogens from MCS when some are matched and others aren't.
    
    Args:
        mol: RDKit molecule
        hits: Tuple of two tuples containing matched atom indices for mol1 and mol2
        mol_index: Integer (0 or 1) indicating which molecule is being processed
    
    Returns:
        Updated hits tuple with symmetrically equivalent hydrogens removed from both molecules
    """
    # Convert tuples to lists so we can modify them
    hits_this_mol = list(hits[mol_index])
    hits_other_mol = list(hits[1 - mol_index])
    
    # Use a set to store H atoms that should be removed (ensures uniqueness)
    h_atoms_to_remove = set()
    
    # Check each atom in the molecule
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() > 1:  # Heavy atom (not H)
            # Get all bonded H atoms
            bonded_h_atoms = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 1]
            
            if len(bonded_h_atoms) >= 2:  # At least two H atoms
                # Check if some H atoms are in MCS and some aren't
                h_indices = [h.GetIdx() for h in bonded_h_atoms] + [atom.GetIdx()]
                h_in_mcs = [idx for idx in h_indices if idx in hits_this_mol]
                h_not_in_mcs = [idx for idx in h_indices if idx not in hits_this_mol]
                
                if h_in_mcs and h_not_in_mcs:  # Some H in MCS, some not
                    h_atoms_to_remove.update(h_indices)  # Add all H atoms to set
    
    # Remove marked H atoms from both molecules' MCS
    for h_idx in h_atoms_to_remove:
        if h_idx in hits_this_mol:
            pos = hits_this_mol.index(h_idx)
            other_idx = hits_other_mol[pos]
            hits_this_mol.remove(h_idx)
            hits_other_mol.remove(other_idx)
    
    # Convert lists back to tuples
    if mol_index == 0:
        return (tuple(hits_this_mol), tuple(hits_other_mol))
    else:
        return (tuple(hits_other_mol), tuple(hits_this_mol))



def exclude_substituents_from_mcs(mol, substituent_info, hits, patterns, mol_index):
    """
    Exclude specific substituents from the MCS list based on given patterns, #keeping anchor atoms.
    
    Args:
        mol: RDKit molecule
        substituent_info: Dictionary from collect_substituent_info containing substituent analysis results
        hits: Tuple of two tuples containing matched atom indices for mol1 and mol2
        patterns: List of patterns to check for exclusion (e.g., [{'anchor': 'C', 'neighbors': ['H', 'H', 'H']}])
        mol_index: Integer (0 or 1) indicating which molecule is being processed
    
    Returns:
        Updated hits tuple with specified substituents excluded from both molecules
    """
    # Convert tuples to lists so we can modify them
    hits_this_mol = list(hits[mol_index])
    hits_other_mol = list(hits[1 - mol_index])
    
    for site_label, data in substituent_info.items():
        anchor_idx = data['anchor_and_neighbors'][0]
        anchor_atom = mol.GetAtomWithIdx(anchor_idx)
        anchor_symbol = anchor_atom.GetSymbol()
        
        # Get all neighbors using RDKit
        neighbors = list(anchor_atom.GetNeighbors())
        neighbor_symbols = [n.GetSymbol() for n in neighbors]
        
        for pattern in patterns:
            if anchor_symbol == pattern['anchor']:
                # Count required specific atoms (excluding '?')
                required_specific = {}
                num_wildcards = 0
                for n in pattern['neighbors']:
                    if n == '?':
                        num_wildcards += 1
                    else:
                        required_specific[n] = required_specific.get(n, 0) + 1
                
                # Count actual occurrences
                actual_counts = {}
                for n in neighbor_symbols:
                    actual_counts[n] = actual_counts.get(n, 0) + 1
                
                # Check if we have enough of each specifically required atom type
                if all(actual_counts.get(n, 0) >= count for n, count in required_specific.items()):
                    # Get indices of matching neighbors
                    matching_neighbor_indices = []
                    remaining_specific = required_specific.copy()
                    
                    # First match specific atoms
                    for neighbor in neighbors:
                        symbol = neighbor.GetSymbol()
                        if symbol in remaining_specific and remaining_specific[symbol] > 0:
                            matching_neighbor_indices.append(neighbor.GetIdx())
                            remaining_specific[symbol] -= 1
                    
                    # Then fill remaining slots with any atoms for wildcards
                    remaining_neighbors = [n for n in neighbors if n.GetIdx() not in matching_neighbor_indices]
                    wildcards_needed = min(num_wildcards, len(remaining_neighbors))
                    
                    if wildcards_needed > 0:
                        for neighbor in remaining_neighbors[:wildcards_needed]:
                            matching_neighbor_indices.append(neighbor.GetIdx())
                    
                    # If we matched all required atoms (specific + wildcards)
                    if len(matching_neighbor_indices) == len(pattern['neighbors']):
#                       # Only exclude matching neighbors from MCS in both molecules (not the anchor)
#                       for idx in matching_neighbor_indices:                        
                        # Exclude anchor and matching neighbors from MCS in both molecules
                        indices_to_exclude = [anchor_idx] + matching_neighbor_indices
                        for idx in indices_to_exclude:
                            if idx in hits_this_mol:
                                pos = hits_this_mol.index(idx)
                                other_idx = hits_other_mol[pos]
                                hits_this_mol.remove(idx)
                                hits_other_mol.remove(other_idx)
                        break  # Stop checking other patterns for this substituent
    
    # Convert lists back to tuples
    if mol_index == 0:
        return (tuple(hits_this_mol), tuple(hits_other_mol))
    else:
        return (tuple(hits_other_mol), tuple(hits_this_mol))


def process_substituent_coord_xyz_zmat(m2sub_xyz_zmat, mol1, hits):
    """
    Process substituent coordinates by finding corresponding anchor atoms and generating new coordinates.
    
    Args:
        m2sub_xyz_zmat (dict): Dictionary containing substituent information and z-matrix
        mol1 (rdkit.Chem.Mol): Reference molecule
        hits (tuple): Tuple of two lists containing matching atom indices (mol1_indices, mol2_indices)
        
    Returns:
        dict: Updated m2sub_xyz_zmat dictionary with new coordinates and reference information
    """
    # Find indices in mol1 of anchor atoms in mol2
    mol2_first3_indices = m2sub_xyz_zmat['indices'][:3]
    mol1_first3_indices = []
    
    for idx in mol2_first3_indices:
        try:
            position = hits[1].index(idx)
            mol1_idx = hits[0][position]
            mol1_first3_indices.append(mol1_idx)
        except ValueError:
            print(f"Warning: Index {idx} from mol2 not found in MCS matches")
            mol1_first3_indices.append(None)
            
    m2sub_xyz_zmat['indices_refmol'] = mol1_first3_indices

    # Get the conformer for 3D coordinates of mol1
    conf = mol1.GetConformer()
    
    # Extract XYZ coordinates for the first 3 indices in 'indices_refmol'
    refmol_anchor_coords = []
    for idx in mol1_first3_indices:
        if idx is not None:
            pos = conf.GetAtomPosition(idx)
            refmol_anchor_coords.append([pos.x, pos.y, pos.z])
        else:
            refmol_anchor_coords.append([None, None, None])
            
    m2sub_xyz_zmat['refmol_anchor_coords'] = refmol_anchor_coords

    # Generate new XYZ coordinates with seed atoms and z-matrix
    atom_names, new_xyz = generate_xyz_from_zmat_with_seeds(
        m2sub_xyz_zmat['zmat_file'], 
        m2sub_xyz_zmat['refmol_anchor_coords']
    )
    assert atom_names == m2sub_xyz_zmat['elements']
    
    # Clean up zmat file
    os.remove(m2sub_xyz_zmat['zmat_file'])
    m2sub_xyz_zmat['zmat_file'] = None
    m2sub_xyz_zmat['coordinates_new'] = new_xyz
    
    return m2sub_xyz_zmat



def generate_zmatrix_file_engine2(atom_info):
    """
    using xyz_to_zmat.py to generate zmatrix file

    Args:
        atom_info: dictionary from collect_atom_info containing:
            - 'elements': list of atom elements
            - 'coordinates': numpy array of XYZ coordinates


    """

    # Create a temporary file
    temp_dir = tempfile.gettempdir()
    temp_xyz_file = tempfile.NamedTemporaryFile(
        prefix='xyz_',
        suffix='.txt',
        mode='w',
        delete=False,
        dir=temp_dir
    )

    temp_zmat_file = tempfile.NamedTemporaryFile(
        prefix='zmat_',
        suffix='.txt',
        mode='w',
        delete=False,
        dir=temp_dir
    )


    coords = atom_info['coordinates']
    elements = atom_info['elements']

    # Generate XYZ file content
    xyz_content = str(len(elements)) + '\n'
    xyz_content += 'title\n'
    for i in range(len(elements)):
        xyz_content += f"{elements[i]}  {coords[i][0]:.6f}  {coords[i][1]:.6f}  {coords[i][2]:.6f}\n"

    # Create temporary XYZ file
    with open(temp_xyz_file.name, 'w') as f:
        f.write(xyz_content)

    # Use xyz_to_zmat to convert to z-matrix
    try:
        xyz_to_zmat(temp_xyz_file.name, temp_zmat_file.name)
        atom_info['zmat_file'] = temp_zmat_file.name
    except Exception as e:
        print(f"Error converting XYZ to Z-matrix: {e}")
        return  # Return without modifying atom_info if conversion fails

 

    return 


# file IO

def read_sdf_coordinates(sdf_file):
    """
    Read only XYZ coordinates from an SDF file.
    
    Args:
        sdf_file (str): Path to SDF file
        
    Returns:
        numpy.ndarray: XYZ coordinates array of shape (n_atoms, 3)
    """
    coordinates = []
    
    with open(sdf_file, 'r') as f:
        lines = f.readlines()
        
        # Get number of atoms from line 4
        n_atoms = int(lines[3].split()[0])
        
        # Read coordinates from atom block
        for i in range(n_atoms):
            line = lines[i + 4]
            coords = [float(x) for x in line.split()[:3]]
            coordinates.append(coords)
    
    return np.array(coordinates)

def modify_sdf_coordinates(sdf_file, new_coordinates, output_file=None):
    """
    Modify XYZ coordinates in an SDF file while preserving all other information.
    
    Args:
        sdf_file (str): Path to input SDF file
        new_coordinates (numpy.ndarray): New XYZ coordinates array of shape (n_atoms, 3)
        output_file (str, optional): Path to output SDF file. If None, will modify input file name
        
    Returns:
        str: Path to the new SDF file
    """
    # Read the original SDF file
    with open(sdf_file, 'r') as f:
        lines = f.readlines()
    
    # Get number of atoms from line 4 (first 3 characters)
    n_atoms = int(lines[3].split()[0])
    
    # Verify coordinate array matches number of atoms
    if len(new_coordinates) != n_atoms:
        raise ValueError(f"Number of coordinates ({len(new_coordinates)}) does not match "
                       f"number of atoms in SDF file ({n_atoms})")
    
    # Create new lines with modified coordinates
    new_lines = []
    atom_section = False
    atom_count = 0
    
    for i, line in enumerate(lines):
        if i < 4:  # Header lines
            new_lines.append(line)
        elif atom_count < n_atoms:  # Atom block
            # Split the line to preserve all other information
            parts = line.split()
            if len(parts) >= 4:  # Valid atom line
                # Format coordinates with same precision as original
                new_line = "{:>10.4f}{:>10.4f}{:>10.4f}".format(
                    new_coordinates[atom_count][0],
                    new_coordinates[atom_count][1],
                    new_coordinates[atom_count][2]
                )
                # Add remaining columns exactly as they were
                remaining = line[30:].rstrip()
                new_line = new_line + remaining + '\n'
                new_lines.append(new_line)
                atom_count += 1
            else:
                new_lines.append(line)
        else:  # Bond block and other data
            new_lines.append(line)
    
    # Determine output file name if not provided
    if output_file is None:
        base, ext = os.path.splitext(sdf_file)
        output_file = f"{base}_modified{ext}"
    
    # Write modified SDF file
    with open(output_file, 'w') as f:
        f.writelines(new_lines)
    
    return output_file





# below are very basic functions 

def setZeroRemoveCoreExtend(AdMat,CoreIndex):
    '''
    This function removes the connection between selected nodes and all other nodes. 
    '''
    natom=len(AdMat)
    ncore=len(CoreIndex)
    
    for i in range(ncore):
        for j in range(natom):
            x=CoreIndex[i]
            y=j
            AdMat[x][y]=0
            AdMat[y][x]=0
            #print(x,y)
    return         

def rmMCSfromConnComp(mcs_index,conn_comp_label):
    
    natom = len(conn_comp_label)
    
    for i in range(len(mcs_index)):
        conn_comp_label[mcs_index[i]]=-1 # set label to -1 for MCS atoms         
    
    return conn_comp_label


def find_most_frequent_elements(list1):
  # Create a dictionary to store the frequency of each element in the list
  freq_dict = {}
  
  # Iterate over the list and update the frequency of each element in the dictionary
  for element in list1:
    if element in freq_dict:
      freq_dict[element] += 1
    else:
      freq_dict[element] = 1

  # Find the maximum frequency in the dictionary
  max_freq = max(freq_dict.values())

  # Iterate over the dictionary and find the elements with the maximum frequency
  most_frequent_elements = []
  for element, freq in freq_dict.items():
    if freq == max_freq:
      most_frequent_elements.append(element)

  # Store the indices of the elements with the maximum frequency in a list
  indices = []
  for i, element in enumerate(list1):
    if element in most_frequent_elements:
      indices.append(i)

  return indices


def get_sublist(list1, indices):
  # Create a new list to store the sublist
  sublist = []

  # Iterate over the list of indices
  for index in indices:
    # Append the corresponding element from the original list to the new list
    sublist.append(list1[index])

  return sublist    

def remove_hydrogens_and_get_mapping_for_sdf(sdf_file):
    """
    Remove hydrogens from an SDF file and track atom mapping between original and new structure.
    
    Args:
        sdf_file (str): Path to input SDF file
        
    Returns:
        tuple: (new_sdf_contents, mapping dictionary {new_idx: old_idx}[0-based])
    """
    with open(sdf_file, 'r') as f:
        lines = f.readlines()
    
    # Get number of atoms and bonds from line 4
    counts_line = lines[3]
    n_atoms = int(counts_line.split()[0])
    n_bonds = int(counts_line.split()[1])
    
    # Process atom block
    atom_lines = lines[4:4+n_atoms]
    bond_lines = lines[4+n_atoms:4+n_atoms+n_bonds]
    remaining_lines = lines[4+n_atoms+n_bonds:]  # Store remaining lines for processing M CHG
    
    # Create mapping and filter out hydrogens
    mapping = {}
    reverse_mapping = {}  # For updating M CHG lines
    new_atom_lines = []
    new_idx = 0
    h_indices = []
    
    for idx, line in enumerate(atom_lines):
        fields = line.split()
        atom_symbol = fields[3]
        
        if atom_symbol != 'H':
            mapping[new_idx] = idx
            reverse_mapping[idx] = new_idx  # Add reverse mapping
            new_atom_lines.append(line)
            new_idx += 1
        else:
            h_indices.append(idx)
    
    # Update bond block
    new_bond_lines = []
    for bond_line in bond_lines:
        fields = bond_line.split()
        atom1_idx = int(fields[0]) - 1
        atom2_idx = int(fields[1]) - 1
        
        if atom1_idx in h_indices or atom2_idx in h_indices:
            continue
            
        new_atom1_idx = reverse_mapping[atom1_idx] + 1  # Convert to 1-based indexing
        new_atom2_idx = reverse_mapping[atom2_idx] + 1
        
        new_bond_line = f"{new_atom1_idx:3d}{new_atom2_idx:3d}{fields[2]:3s}{fields[3]:3s}\n"
        new_bond_lines.append(new_bond_line)
    
    # Process M CHG lines and other remaining content
    new_remaining_lines = []
    i = 0
    while i < len(remaining_lines):
        line = remaining_lines[i]
        if line.startswith('M  CHG'):
            fields = line.split()
            n_charges = int(fields[2])  # Number of charge entries
            new_charges = []
            
            # Process each charge entry
            for j in range(n_charges):
                atom_idx = int(fields[3 + j*2]) - 1  # Convert to 0-based
                charge = fields[4 + j*2]
                
                # Only keep charges for non-hydrogen atoms
                if atom_idx not in h_indices:
                    new_atom_idx = reverse_mapping[atom_idx] + 1  # Convert back to 1-based
                    new_charges.extend([str(new_atom_idx), charge])
            
            if new_charges:  # Only add M CHG line if there are charges to record
                new_chg_line = f"M  CHG  {len(new_charges)//2}"
                for val in new_charges:
                    new_chg_line += f" {val:>3s}"
                new_chg_line += "\n"
                new_remaining_lines.append(new_chg_line)
        else:
            new_remaining_lines.append(line)
        i += 1
    
    # Create new SDF content
    new_lines = lines[:3]  # Keep header
    new_counts_line = f"{len(new_atom_lines):3d}{len(new_bond_lines):3d}  0  0  0  0  0  0  0  0999 V2000\n"
    new_lines.append(new_counts_line)
    new_lines.extend(new_atom_lines)
    new_lines.extend(new_bond_lines)
    new_lines.extend(new_remaining_lines)
    
    return ''.join(new_lines), mapping

def write_processed_sdf(sdf_contents, output_file):
    """
    Write the processed SDF contents to a file.
    
    Args:
        sdf_contents (str): SDF file contents as string
        output_file (str): Path to output SDF file
    """
    with open(output_file, 'w') as f:
        f.write(sdf_contents)

