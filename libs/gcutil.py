# This file is adapted from the original gcutil.py file provided by Robert A Shaw
# https://github.com/robashaw/geomConvert/blob/master/gcutil.py

# UTILITY TO CONVERT BETWEEN XYZ AND Z-MATRIX GEOMETRIES
# Copyright 2017 Robert A Shaw
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the Software
# is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#                                                   
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
# OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
# Utilities for gc.py



import numpy as np
from scipy.spatial.distance import cdist
import tempfile
import os
import sys


def replace_vars(vlist, variables):
    """ Replaces a list of variable names (vlist) with their values
        from a dictionary (variables).
    """
    for i, v in enumerate(vlist):
        if v in variables:
            vlist[i] = variables[v]
        else:
            try:
                # assume the "variable" is a number
                vlist[i] = float(v)
            except:
                print("Problem with entry " + str(v))

def readxyz(filename):
    """ Reads in a .xyz file in the standard format,
        returning xyz coordinates as a numpy array
        and a list of atom names.
    """
    xyzf = open(filename, 'r')
    xyzarr = np.zeros([1, 3])
    atomnames = []
    if not xyzf.closed:
        # Read the first line to get the number of particles
        npart = int(xyzf.readline())
        # and next for title card
        title = xyzf.readline()

        # Make an N x 3 matrix of coordinates
        xyzarr = np.zeros([npart, 3])
        i = 0
        for line in xyzf:
            words = line.split()
            if (len(words) > 3):
                atomnames.append(words[0])
                xyzarr[i][0] = float(words[1])
                xyzarr[i][1] = float(words[2])
                xyzarr[i][2] = float(words[3])
                i = i + 1
    return (xyzarr, atomnames)

def readzmat(filename):
    """ Reads in a z-matrix in standard format,
        returning a list of atoms and coordinates.
    """
    zmatf = open(filename, 'r')
    atomnames = []
    rconnect = []  # bond connectivity
    rlist = []     # list of bond length values
    aconnect = []  # angle connectivity
    alist = []     # list of bond angle values
    dconnect = []  # dihedral connectivity
    dlist = []     # list of dihedral values
    variables = {} # dictionary of named variables
    
    if not zmatf.closed:
        for line in zmatf:
            words = line.split()
            eqwords = line.split('=')
            
            if len(eqwords) > 1:
                # named variable found 
                varname = str(eqwords[0]).strip()
                try:
                    varval  = float(eqwords[1])
                    variables[varname] = varval
                except:
                    print("Invalid variable definition: " + line)
            
            else:
                # no variable, just a number
                # valid line has form
                # atomname index1 bond_length index2 bond_angle index3 dihedral
                if len(words) > 0:
                    atomnames.append(words[0])
                if len(words) > 1:
                    rconnect.append(int(words[1]))
                if len(words) > 2:
                    rlist.append(words[2])
                if len(words) > 3:
                    aconnect.append(int(words[3]))
                if len(words) > 4:
                    alist.append(words[4])
                if len(words) > 5:
                    dconnect.append(int(words[5]))
                if len(words) > 6:
                    dlist.append(words[6])
    
    # replace named variables with their values
    replace_vars(rlist, variables)
    replace_vars(alist, variables)
    replace_vars(dlist, variables)
    
    return (atomnames, rconnect, rlist, aconnect, alist, dconnect, dlist) 

def distance_matrix(xyzarr):
    """Returns the pairwise distance matrix between atom
       from a set of xyz coordinates 
    """
    return cdist(xyzarr, xyzarr)

def angle(xyzarr, i, j, k):
    """Return the bond angle in degrees between three atoms 
       with indices i, j, k given a set of xyz coordinates.
       atom j is the central atom
    """
    rij = xyzarr[i] - xyzarr[j]
    rkj = xyzarr[k] - xyzarr[j]
    cos_theta = np.dot(rij, rkj)
    sin_theta = np.linalg.norm(np.cross(rij, rkj))
    theta = np.arctan2(sin_theta, cos_theta)
    theta = 180.0 * theta / np.pi 
    return theta

def dihedral(xyzarr, i, j, k, l):
    """Return the dihedral angle in degrees between four atoms 
       with indices i, j, k, l given a set of xyz coordinates.
       connectivity is i->j->k->l
    """
    rji = xyzarr[j] - xyzarr[i]
    rkj = xyzarr[k] - xyzarr[j]
    rlk = xyzarr[l] - xyzarr[k]
    v1 = np.cross(rji, rkj)
    v1 = v1 / np.linalg.norm(v1)
    v2 = np.cross(rlk, rkj)
    v2 = v2 / np.linalg.norm(v2)
    m1 = np.cross(v1, rkj) / np.linalg.norm(rkj)
    x = np.dot(v1, v2)
    y = np.dot(m1, v2)
    chi = np.arctan2(y, x)
    chi = -180.0 - 180.0 * chi / np.pi
    if (chi < -180.0):
        chi = chi + 360.0
    return chi

def write_zmat(xyzarr, distmat, atomnames, rvar=False, avar=False, dvar=False):
    """Prints a z-matrix from xyz coordinates, distances, and atomnames,
       optionally with the coordinate values replaced with variables.
    """
    npart, ncoord = xyzarr.shape
    rlist = [] # list of bond lengths
    alist = [] # list of bond angles (degrees)
    dlist = [] # list of dihedral angles (degrees)
    if npart > 0:
        # Write the first atom
        print(atomnames[0])
        
        if npart > 1:
            # and the second, with distance from first
            n = atomnames[1]
            rlist.append(distmat[0][1])
            if (rvar):
                r = 'R1'
            else:
                r = '{:>11.5f}'.format(rlist[0])
            print('{:<3s} {:>4d}  {:11s}'.format(n, 1, r))
            
            if npart > 2:
                n = atomnames[2]
                
                rlist.append(distmat[0][2])
                if (rvar):
                    r = 'R2'
                else:
                    r = '{:>11.5f}'.format(rlist[1])
                
                alist.append(angle(xyzarr, 2, 0, 1))
                if (avar):
                    t = 'A1'
                else:
                    t = '{:>11.5f}'.format(alist[0])

                print('{:<3s} {:>4d}  {:11s} {:>4d}  {:11s}'.format(n, 1, r, 2, t))
                
                if npart > 3:
                    for i in range(3, npart):
                        n = atomnames[i]

                        rlist.append(distmat[i-3][i])
                        if (rvar):
                            r = 'R{:<4d}'.format(i)
                        else:
                            r = '{:>11.5f}'.format(rlist[i-1])

                        alist.append(angle(xyzarr, i, i-3, i-2))
                        if (avar):
                            t = 'A{:<4d}'.format(i-1)
                        else:
                            t = '{:>11.5f}'.format(alist[i-2])
                        
                        dlist.append(dihedral(xyzarr, i, i-3, i-2, i-1))
                        if (dvar):
                            d = 'D{:<4d}'.format(i-2)
                        else:
                            d = '{:>11.5f}'.format(dlist[i-3])
                        print('{:3s} {:>4d}  {:11s} {:>4d}  {:11s} {:>4d}  {:11s}'.format(n, i-2, r, i-1, t, i, d))
    if (rvar):
        print(" ")
        for i in range(npart-1):
            print('R{:<4d} = {:>11.5f}'.format(i+1, rlist[i]))
    if (avar):
        print(" ")
        for i in range(npart-2):
            print('A{:<4d} = {:>11.5f}'.format(i+1, alist[i]))
    if (dvar):
        print(" ")
        for i in range(npart-3):
            print('D{:<4d} = {:>11.5f}'.format(i+1, dlist[i]))

def generate_zmatrix_file(atom_info, rvar=False, avar=False, dvar=False, zmat_gen_method=1):
    """
    Generate Z-matrix from atom elements and coordinates and write to a temporary file.
    
    Args:
        atom_info: dictionary from collect_atom_info containing:
            - 'elements': list of atom elements
            - 'coordinates': numpy array of XYZ coordinates
        rvar, avar, dvar: booleans to control variable naming in Z-matrix
        zmat_gen_method: int, method to generate z-matrix (1: use our own zmat generation, 2: use geomConvert's method)
            
    Returns:
        str: Path to temporary file containing the Z-matrix
    """
    # Create a temporary file
    temp_dir = tempfile.gettempdir()
    temp_file = tempfile.NamedTemporaryFile(
        prefix='zmatrix_',
        suffix='.txt',
        mode='w',
        delete=False,
        dir=temp_dir
    )
    
    # Redirect stdout to the temporary file
    import sys
    original_stdout = sys.stdout
    sys.stdout = temp_file
    
    # Get distance matrix between all atoms
    coords = atom_info['coordinates']
    elements = atom_info['elements']
    distmat = distance_matrix(coords)
    
    # Generate the Z-matrix
    if zmat_gen_method == 2: #gcutil's write_zmat function   
        write_zmat(coords, distmat, elements, rvar=rvar, avar=avar, dvar=dvar)
    elif zmat_gen_method == 1: # our own zmat generation
        write_zmat_by_connectivity(coords, distmat, elements, rvar=rvar, avar=avar, dvar=dvar)
    else:
        raise ValueError('Invalid zmat generation method')
    
    # Restore stdout and close file
    sys.stdout = original_stdout
    temp_file.close()
    
    atom_info['zmat_file'] = temp_file.name 
    
    return


def write_xyz(atomnames, rconnect, rlist, aconnect, alist, dconnect, dlist, first3_coords=None):
    """Prints out an xyz file from a decomposed z-matrix
    Args:
        ...existing args...
        first3_coords: numpy array of shape (3,3) containing XYZ coordinates for first 3 atoms
    """
    npart = len(atomnames)
    print(npart)
    print('INSERT TITLE CARD HERE')
    
    xyzarr = np.zeros([npart, 3])
    
    if first3_coords is not None:
        # Use provided coordinates for first 3 atoms
        xyzarr[0:3] = first3_coords
    else:
        # Original code for first 3 atoms
        if (npart > 1):
            xyzarr[1] = [rlist[0], 0.0, 0.0]
        if (npart > 2):
            i = rconnect[1] - 1
            j = aconnect[0] - 1
            r = rlist[1]
            theta = alist[0] * np.pi / 180.0
            x = r * np.cos(theta)
            y = r * np.sin(theta)
            a_i = xyzarr[i]
            b_ij = xyzarr[j] - xyzarr[i]
            if (b_ij[0] < 0):
                x = a_i[0] - x
                y = a_i[1] - y
            else:
                x = a_i[0] + x
                y = a_i[1] + y
            xyzarr[2] = [x, y, 0.0]

    for n in range(3, npart):
        # back-compute the xyz coordinates
        # from the positions of the last three atoms
        r = rlist[n-1]
        theta = alist[n-2] * np.pi / 180.0
        phi = dlist[n-3] * np.pi / 180.0
        
        sinTheta = np.sin(theta)
        cosTheta = np.cos(theta)
        sinPhi = np.sin(phi)
        cosPhi = np.cos(phi)

        x = r * cosTheta
        y = r * cosPhi * sinTheta
        z = r * sinPhi * sinTheta
        
        i = rconnect[n-1] - 1
        j = aconnect[n-2] - 1
        k = dconnect[n-3] - 1
        a = xyzarr[k]
        b = xyzarr[j]
        c = xyzarr[i]
        
        ab = b - a
        bc = c - b
        bc = bc / np.linalg.norm(bc)
        nv = np.cross(ab, bc)
        nv = nv / np.linalg.norm(nv)
        ncbc = np.cross(nv, bc)
        
        new_x = c[0] - bc[0] * x + ncbc[0] * y + nv[0] * z
        new_y = c[1] - bc[1] * x + ncbc[1] * y + nv[1] * z
        new_z = c[2] - bc[2] * x + ncbc[2] * y + nv[2] * z
        xyzarr[n] = [new_x, new_y, new_z]
            
    # print results
    for i in range(npart):
        print('{:<4s}\t{:>11.5f}\t{:>11.5f}\t{:>11.5f}'.format(atomnames[i], xyzarr[i][0], xyzarr[i][1], xyzarr[i][2]))


def generate_xyz_from_zmat_with_seeds(zmat_file, seed_coords):
    """
    Generate XYZ coordinates using a Z-matrix file and seed coordinates for the first three atoms.
    
    Args:
        zmat_file (str): Path to Z-matrix file
        seed_coords (list): List of [x,y,z] coordinates for the first three atoms
        
    Returns:
        tuple: (atomnames, xyz_coords) where xyz_coords is a numpy array of shape (n_atoms, 3)
    """
    # Read the Z-matrix file
    atomnames, rconnect, rlist, aconnect, alist, dconnect, dlist = readzmat(zmat_file)
    
    # Convert seed coordinates to numpy array
    seed_coords = np.array(seed_coords)
    
    # Create temporary file for output
    temp_dir = tempfile.gettempdir()
    temp_file = tempfile.NamedTemporaryFile(
        prefix='xyz_',
        suffix='.txt',
        mode='w',
        delete=False,
        dir=temp_dir
    )
    
    # Redirect stdout to capture write_xyz output
    original_stdout = sys.stdout
    sys.stdout = temp_file
    
    # Generate XYZ coordinates using the seed coordinates
    write_xyz(atomnames, rconnect, rlist, aconnect, alist, dconnect, dlist, first3_coords=seed_coords)
    
    # Restore stdout and close file
    sys.stdout = original_stdout
    temp_file.close()
    
    # Read the generated XYZ file
    xyz_coords, atom_names = readxyz(temp_file.name)
    
    # Clean up temporary file
    os.remove(temp_file.name)
    
    return atom_names, xyz_coords

def write_zmat_by_connectivity(xyzarr, distmat, atomnames, rvar=False, avar=False, dvar=False):
    """Prints a z-matrix from xyz coordinates using actual molecular connectivity.
    For each atom after the third one:
    1. Bond partner (j): closest atom among previous atoms
    2. Angle partner (k): bond partner of j
    3. Dihedral partner (l): bond partner of k
    
    Rules:
    - Follow bond partner chain unless atom 1 is encountered
    - If atom 1 appears as any reference, other references must be from atoms 2 and 3
    - All three reference atoms must be different
    """
    npart, ncoord = xyzarr.shape
    rlist = []  # list of bond lengths
    alist = []  # list of bond angles (degrees)
    dlist = []  # list of dihedral angles (degrees)
    
    # Store connectivity information for each atom
    bond_partners = {}  # key: atom index, value: its bond partner's index
    
    if npart > 0:
        # First atom
        print(atomnames[0])
        
        if npart > 1:
            # Second atom - always connected to first
            n = atomnames[1]
            rlist.append(distmat[0][1])
            bond_partners[1] = 0
            if (rvar):
                r = 'R1'
            else:
                r = '{:>11.5f}'.format(rlist[0])
            print('{:<3s} {:>4d}  {:11s}'.format(n, 1, r))
            
            if npart > 2:
                # Third atom - connected to first atom by default
                n = atomnames[2]
                bond_partners[2] = 0
                
                rlist.append(distmat[0][2])
                if (rvar):
                    r = 'R2'
                else:
                    r = '{:>11.5f}'.format(rlist[1])
                
                alist.append(angle(xyzarr, 2, 0, 1))
                if (avar):
                    t = 'A1'
                else:
                    t = '{:>11.5f}'.format(alist[0])
                
                print('{:<3s} {:>4d}  {:11s} {:>4d}  {:11s}'.format(n, 1, r, 2, t))
                
                if npart > 3:
                    for i in range(3, npart):
                        n = atomnames[i]
                        
                        # Find closest previous atom (bond partner j)
                        prev_dists = distmat[0:i, i]
                        j = np.argmin(prev_dists)
                        bond_partners[i] = j
                        
                        # Get k by following bond partner chain
                        k = bond_partners.get(j, 0)
                        
                        # Get l by following bond partner chain
                        l = bond_partners.get(k, 0)
                        
                        # Now adjust references if atom 1 is involved
                        if j == 0:
                            # If directly bonded to atom 1
                            k = 1
                            l = 2
                        elif k == 0:
                            # If second reference is atom 1
                            if j == 1:
                                l = 2
                            elif j == 2:
                                l = 1
                            else:
                                l = 2
                        
                        # Ensure all references are different
                        if len(set([j, k, l])) < 3:
                            # If we have duplicates, adjust l
                            for possible_l in [1, 2, 3]:
                                if possible_l not in [j, k]:
                                    l = possible_l
                                    break
                        
                        # Calculate values and print
                        rlist.append(distmat[j][i])
                        if (rvar):
                            r = 'R{:<4d}'.format(i)
                        else:
                            r = '{:>11.5f}'.format(rlist[i-1])
                        
                        alist.append(angle(xyzarr, i, j, k))
                        if (avar):
                            t = 'A{:<4d}'.format(i-1)
                        else:
                            t = '{:>11.5f}'.format(alist[i-2])
                        
                        dlist.append(dihedral(xyzarr, i, j, k, l))
                        if (dvar):
                            d = 'D{:<4d}'.format(i-2)
                        else:
                            d = '{:>11.5f}'.format(dlist[i-3])
                        
                        # Print atom line with 1-based indices
                        print('{:3s} {:>4d}  {:11s} {:>4d}  {:11s} {:>4d}  {:11s}'.format(
                            n, j+1, r, k+1, t, l+1, d))
    
    # Print variable definitions if requested
    if (rvar):
        print(" ")
        for i in range(npart-1):
            print('R{:<4d} = {:>11.5f}'.format(i+1, rlist[i]))
    if (avar):
        print(" ")
        for i in range(npart-2):
            print('A{:<4d} = {:>11.5f}'.format(i+1, alist[i]))
    if (dvar):
        print(" ")
        for i in range(npart-3):
            print('D{:<4d} = {:>11.5f}'.format(i+1, dlist[i]))
