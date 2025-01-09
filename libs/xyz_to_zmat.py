"""
XYZ-to-ZMAT.py from https://github.com/tjmustard/XYZtoZMAT/

This script converts an XYZ file to a Z-matrix format.
"""

import os
import math
from decimal import Decimal
import sys
import numpy as np
import argparse
import getopt

def xyz_to_zmat(input_file, output_file, print_xyz=False):
	"""
	Converts an XYZ file to a Z-matrix format.

	Args:
		input_file: Path to the input XYZ file.
		output_file: Path to the output Z-matrix file.
		print_xyz: Boolean, whether to print the XYZ coordinates. Defaults to False.
	"""

	### --- Arguments --- ###
	# Removed global variables, using function arguments instead
	# program = "XYZ-to-ZMAT.py"
	# ifile = ''
	# ofile = ''
	# printfile = 0

	# ifile = input_file
	# ofile = output_file
	# printfile = print_xyz

	### --- Distance function (getdistance) --- ###
	def getdistance(at1, at2, lol):
		try:
			atom1 = np.array([float(lol[at1][1]), float(lol[at1][2]), float(lol[at1][3])])
			atom2 = np.array([float(lol[at2][1]), float(lol[at2][2]), float(lol[at2][3])])
			vector1 = atom2 - atom1
			dist = np.linalg.norm(vector1)
			return dist
		except:
			return 0

	### --- Angle function --- ###
	def getangle(at1, at2, at3, lol):
		try:
			atom1 = np.array([float(lol[at1][1]), float(lol[at1][2]), float(lol[at1][3])])
			atom2 = np.array([float(lol[at2][1]), float(lol[at2][2]), float(lol[at2][3])])
			atom3 = np.array([float(lol[at3][1]), float(lol[at3][2]), float(lol[at3][3])])
			vector1 = atom1 - atom2
			vector2 = atom3 - atom2
			angle = np.arccos(np.dot(vector1, vector2) / (np.linalg.norm(vector1) * np.linalg.norm(vector2)))
			return math.degrees(angle)
		except:
			return 0.000000

	### --- Dihedral angle function --- ###
	def getdihedral(at1, at2, at3, at4, lol):
		try:
			atom1 = np.array([float(lol[at1][1]), float(lol[at1][2]), float(lol[at1][3])])
			atom2 = np.array([float(lol[at2][1]), float(lol[at2][2]), float(lol[at2][3])])
			atom3 = np.array([float(lol[at3][1]), float(lol[at3][2]), float(lol[at3][3])])
			atom4 = np.array([float(lol[at4][1]), float(lol[at4][2]), float(lol[at4][3])])
			vector1 = atom2 - atom1
			vector2 = atom3 - atom2
			plane1 = np.cross(vector1, vector2)
			vector3 = atom2 - atom3
			vector4 = atom4 - atom3
			plane2 = np.cross(vector3, vector4)
			dihedral = np.arccos(-np.dot(plane1, plane2) / (np.linalg.norm(plane1) * np.linalg.norm(plane2)))
			if np.dot(plane1, vector4) > 0:
				return math.degrees(dihedral)
			else:
				return -math.degrees(dihedral)
		except:
			return 0

	### --- Functions to get and give element numbers and names --- ###
	elementList = ["h","he","li","be","b","c","n","o","f","ne","na","mg","al","si","p","s","cl","ar","k","ca","sc","ti","v","cr","mn","fe","co","ni","cu","zn","ga","ge","as","se","br","kr","rb","sr","y","zr","nb","mo","tc","ru","rh","pd","ag","cd","in","sn","sb","te","i","xe","cs","ba","la","ce","pr","nd","pm","sm","eu","gd","tb","dy","ho","er","tm","yb","lu","hf","ta","w","re","os","ir","pt","au","hg","tl","pb","bi","po","at","rn","fr","ra","ac","th","pa","u","np","pu","am","cm","bk","cf","es","fm","md","no","lr","rf","db","sg","bh","hs","mt","ds","rg","cn","uut","fl","uup","lv","uus","uuo"]
	elementNames = ["H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Uut","Fl","Uup","Lv","Uus","Uuo"]
	elementLarge =  ["na","mg","al","si","p","s","cl","ar","k","ca","sc","ti","v","cr","mn","fe","co","ni","cu","zn","ga","ge","as","se","br","kr","rb","sr","y","zr","nb","mo","tc","ru","rh","pd","ag","cd","in","sn","sb","te","i","xe","cs","ba","la","ce","pr","nd","pm","sm","eu","gd","tb","dy","ho","er","tm","yb","lu","hf","ta","w","re","os","ir","pt","au","hg","tl","pb","bi","po","at","rn","fr","ra","ac","th","pa","u","np","pu","am","cm","bk","cf","es","fm","md","no","lr","rf","db","sg","bh","hs","mt","ds","rg","cn","uut","fl","uup","lv","uus","uuo"]
	def getElementNum(at1):
		element = elementList.index(at1.lower())
		return element+1
	def getElementName(at1):
		element = elementNames[at1-1]
		return element  # Return element name

	### --- Open parent file --- ###
	try:
		f = open(input_file)
		ifileList = f.readlines()
		f.close()
	except FileNotFoundError:
		print(f"Error: Input file '{input_file}' not found.")
		return

	### --- Create some variables needed for parsing the input file --- ###
	# ifilelength is used to gauge what line we are on and also allows us to determine the number of lines in the file
	# ifilelol is a list of lists (lol) that will hold all the data of the input file
	# a_list is a temporary list that will be appended to the ifilelol
	################################
	ifileLength = 0
	ifilelol = []
	a_list = []
	#### --- Input/parse the input file into a list of lists called ifilelol --- ###
	for i, var in enumerate(ifileList):
		line = var.split()
		#print line
		if i == 0:
			a_list.append(var)
			ifilelol.append(list(a_list))
			a_list[:] =[]
		if i == 1:
			a_list.append(var.rstrip('\n'))
			ifilelol.append(list(a_list))
			a_list[:] =[]
		if i >= 2:
			a_list.append(getElementNum(line[0]))
			a_list.append(Decimal(line[1]))
			a_list.append(Decimal(line[2]))
			a_list.append(Decimal(line[3]))
			ifilelol.append(list(a_list))
			a_list[:] =[]
		ifileLength = ifileLength + 1

	### --- Get bonds --- ###
	iFileAtomNum = ' '.join(ifilelol[0])
	iFileName = ' '.join(ifilelol[1])  # Get the file name
	covBonds = []
	covHBonds = []
	covTMBonds = []
	nearestNeighbor = []
	neighborStart = [0, 1000000]
	#### --- Generate bond lists --- ###
	#for i in range(0,int(iFileAtomNum)):
	#	nearestNeighbor.append(list(neighborStart))
	#	for j in range(0,int(iFileAtomNum)):
	#		if i !=j:
	#			distij = getdistance(i+2, j+2)
	#			if j > i:
	#				if distij <= 2.25 and ifilelol[i+2][0] != 1 and ifilelol[j+2][0] != 1 and ((getElementName(ifilelol[i+2][0]).lower() not in elementLarge) and (getElementName(ifilelol[j+2][0]).lower() not in elementLarge)):
	#					distList = [i+1, j+1, distij]
	#					covBonds.append(distList)
	#					#print str(i+2) + "\t" + str(j+2) + "\t" + str(distij)
	#				elif distij <= 1.3 and (ifilelol[i+2][0] == 1 or ifilelol[j+2][0] == 1):
	#					distList = [i+1, j+1, distij]
	#					covHBonds.append(distList)
	#				elif distij <= 3 and ((getElementName(ifilelol[i+2][0]).lower() in elementLarge) and (getElementName(ifilelol[j+2][0]).lower() in elementLarge)):
	#					distList = [i+1, j+1, distij]
	#					covTMBonds.append(distList)
	#			if distij < nearestNeighbor[i][1]:
	#				nearestNeighbor[i][0] = j + 1
	#				nearestNeighbor[i][1] = distij
	#### --- Remove hydrogen bonds from bond list --- ###
	#for i in range(0,len(covHBonds)):
	#	if (covHBonds[i][0] != nearestNeighbor[covHBonds[i][0]][0]) and (covHBonds[i][1] != nearestNeighbor[covHBonds[i][0]][0]):
	#		del covHBonds[i]	
	##print "Covalent bonds to Hydrogen:"
	##print covHBonds
	##print "Covalent bonds between \"heavy\" atoms."
	##print covBonds
	##print "Covalent bonds between TM atoms."
	##print covTMBonds
	##print nearestNeighbor

	### --- get the distance, angle and dihedrals for a Z-matrix --- ###
	def getzmat(i):
		line = []
		line.append(getElementName(ifilelol[i+2][0]))
		if i > 0:
			line.append(i)
			dist = getdistance(i+1, i+2, ifilelol)
			line.append(dist)
		if i > 1:
			line.append(i-1)
			angle = getangle(i, i+1, i+2, ifilelol)
			line.append(angle)
		if i > 2:
			line.append(i-2)
			dihedral = getdihedral(i-1, i, i+1, i+2, ifilelol)
			line.append(dihedral)
		line.append(-1)
		line.append(-1)
		line.append(-1)
		line.append(-1)
		line.append(-1)
		line.append(-1)
		line.append("\n")
		return line

	### --- Get the XYZ coordinates from distance, angle and dihedral data --- ###
	def getXYZfromZMAT(lol, at4):
		### Set up the variables to be used for the function
		dist = float(lol[2])
		angle = float(lol[4])
		dihedral = float(lol[6])
		angleRad = math.radians(angle) # * math.pi / 180
		dihedralRad = math.radians(dihedral) # * math.pi / 180
		at1 = int(lol[5])-1
		at2 = int(lol[3]) - 1
		at3 = int(lol[1])-1
		x = 0
		y = 0
		z = 0
		line = []

		### Start to place the atoms in their locations
		if at4 == 0:
			x = 0.00000
			y = 0.00000
			z = 0.00000
		
		elif at4 == 1:  # Handle the second atom
			x = dist
			y = 0.00000
			z = 0.00000
		
		elif at4 == 2:
			a = xyzLOL[at3][1]
			b = dist
			x = a + (dist * math.cos(math.pi - angleRad))
			y = 0.00000
			z = -dist * math.sin(math.pi - angleRad)
		
		elif at4 >= 3:
			####The at4 x,y,z coordinates from spherical coord
			Sx = dist * math.sin(angleRad) * math.cos(dihedralRad)
			Sy = -dist * math.sin(angleRad) * math.sin(dihedralRad) #For some reason I need to have a negative here to get the correct sign in the output..... weird
			Sz = dist * math.cos(angleRad)
			at4L = [Sx, Sy, Sz]
			#print "at4L: " + str(at4L)

			###Finding the angle theta
			#Make the list of lists for the three point (z-axis, origin, and translated atom 2) needed for an angle calculation
			Z32 = [[0, 0, 0, 1],  [0, 0, 0, 0],  [0, xyzLOL[at2][1] - xyzLOL[at3][1], xyzLOL[at2][2] - xyzLOL[at3][2], xyzLOL[at2][3] - xyzLOL[at3][3]]]
			#Get theta using the getangle function
			theta = math.radians(getangle(0, 1, 2, Z32))
			#print "theta: " +str(theta)
			
			###Rodrigues' rotation formula
			#Create the vectprs needed to calculate k
			vector3 = np.array([xyzLOL[at3][1], xyzLOL[at3][2], xyzLOL[at3][3]])
			vector2 = np.array([xyzLOL[at2][1], xyzLOL[at2][2], xyzLOL[at2][3]])
			vector0 = np.array([0, 0, 1])
			#Calculate k for the Rodrigues rotation formula
			k = np.cross((vector2-vector3), vector0)/np.linalg.norm(np.cross((vector2-vector3), vector0))
			#Generate an array for translated 1
			T1 = [(xyzLOL[at1][1]-xyzLOL[at3][1]), (xyzLOL[at1][2]-xyzLOL[at3][2]), (xyzLOL[at1][3]-xyzLOL[at3][3])]
			#Calculate the Rodrigues rotation matrix
			RR23T1 = np.dot(T1, np.cos(theta)) + np.dot(np.cross(k,T1), np.sin(theta)) + np.dot(np.dot(k,(np.dot(k,T1))), (1-np.cos(theta)))
			#Make the list of lists for the four points (x-axis, z-axis, origin, and rotated translated 1) needed for a dihedral calculation
			XZ31 = [[0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 0, 0], [0, RR23T1[0], RR23T1[1], RR23T1[2]]]
			#Get phi using the getdihedral function
			phi = math.radians(getdihedral(0,1,2,3,XZ31))
			#print "k: " + str(k)
			#print "T1: " + str(T1)
			#print "RR23T1: " + str(RR23T1)
			#print "XZ31: " + str(XZ31)
			#print("phi: " + str(phi))

			###Rotation matrix 
			#Create the array for the rotation matrix including dihedral phi
			RM = np.array([[np.cos(phi), np.sin(phi), 0], [-np.sin(phi), np.cos(phi), 0], [0, 0, 1]])
			#Calculate the dot product of the rotation matrix and the coordinates for 4 (from spherical)
			
			RM4 = np.dot(RM, at4L)
			#Calculate the rotated coordinates of the rotated coordinates of atom 4
			RRN23RM4 = np.dot(RM4, np.cos(-theta)) + np.dot(np.cross(k,RM4), np.sin(-theta)) + np.dot(np.dot(k,(np.dot(k,RM4))), (1-np.cos(-theta)))
			#print "RM: " + str(RM)
			#print "RM4: " + str(RM4)
			#print "RRN23RM4: " + str(RRN23RM4)
			#Final coordinates that are rotated, rotated and translated
			x = RRN23RM4[0] + xyzLOL[at3][1]
			y = RRN23RM4[1] + xyzLOL[at3][2]
			z = RRN23RM4[2] + xyzLOL[at3][3]
			#print x
			#print y
			#print z
		#Putting everything into a list to send back
		line.append(lol[0])
		line.append(x)
		line.append(y)
		line.append(z)
		#line.append("\n")
		return line

	### --- Generate a list of lists for the zmatrix of the structure --- ###
	zmat_lol = []
	for i in range(0,int(iFileAtomNum)):
		linetemp = getzmat(i)
		zmat_lol.append(linetemp)
		#a_list[:] =[]

	### --- Output file in Z-matrix format --- ###
	f = open(output_file, 'w+')
	#f.write(str(iFileAtomNum))
	#f.write(str(iFileName) + '\n')
	for i in range(0,int(iFileAtomNum)):
		linetemp = [x for x in zmat_lol[i] if x != -1]
		line = '\t'.join(str(e) for e in linetemp) #, e != -1)
		f.write(line)
	f.close()

	### --- Print out the xyz from the zmatrix --- ###
	if print_xyz == True:
		xyzLOL = []
		print(iFileAtomNum)
		#print ""
		for i in range(0,int(iFileAtomNum)):
			linetemp = getXYZfromZMAT(zmat_lol[i],i)
			xyzLOL.append(linetemp)
			print('   '.join(str(e) for e in linetemp))


	######################################################################
	### END OF SCRIPT	
	######################################################################