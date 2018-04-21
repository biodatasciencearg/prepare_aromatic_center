#!/usr/bin/python2.7 

import sys, os
import openbabel, pybel
from itertools import chain
from os.path import basename
import numpy as np
################################################################
####################auxiliary programs##########################
################################################################

def center_of_mass(Vx, Vy, Vz, mass):
	cgx = np.sum(Vx*mass)/np.sum(mass)
	cgy = np.sum(Vy*mass)/np.sum(mass)
	cgz = np.sum(Vz*mass)/np.sum(mass)
	return (cgx, cgy, cgz)


################################################################
####################pdbqt parser#####################
################################################################
	
class PDBQT():
    def __init__(self, line):
        self._parse_common(line)     # used here (PDB) and in PDBQT
        self._parse_specific(line)   # differs in PDBQT

    def getline(self):
        txt = self._print_common()    # no \n; PDB + PDBQT
        txt += self._print_specific() # differs in PDBQT
        return txt
        
    def _parse_common(self, line):
        """Common to PDB and PDBQT formats"""
        self.keyword     = line      [ 0: 6]     # ATOM or HETATM
        self.serial      = int(line  [ 6:11])    # atom id
        #                            [11:12]
        self.name        = line      [12:16]     # atom name
        self.altLoc      = line      [16:17]     # Alternate location
        self.resName     = line      [17:20]     # Residue name
        #                            [20:21] 
        self.chain       = line      [21:22]     # chain
        self.resNum      = int(line  [22:26])    # Residue number
        self.icode       = line      [26:27]     # ???
        #                            [27:30]
        self.x           = float(line[30:38])    # X
        self.y           = float(line[38:46])    # Y
        self.z           = float(line[46:54])    # Z
        self.occupancy   = float(line[54:60])    # Occupancy
        self.bfact       = float(line[60:66])    # Temperature factor

    def _parse_specific(self, line):
        """ PDBQT characters [68:79] """
        self.charge      = float(line[68:76])   # Charge
        self.atype       = line      [77:79]    # Atom type
        self.atype = self.atype.strip().upper()

    def _print_common(self):
        """ Characters [0:68]"""
        linestr = ''
        linestr += '%6s' % (self.keyword)
        linestr += '%5d' % (self.serial)
        linestr += ' ' 
        linestr += '%4s' % (self.name)
        linestr += '%1s' % (self.altLoc) 
        linestr += '%3s' % (self.resName)
        linestr += ' ' 
        linestr += '%1s' % (self.chain)
        linestr += '%4d' % (self.resNum)
        linestr += '%1s' % (self.icode)
        linestr += ' ' * 3 
        linestr += '%8.3f' % (self.x)
        linestr += '%8.3f' % (self.y)
        linestr += '%8.3f' % (self.z)
        linestr += '%6.2f' % (self.occupancy)
        linestr += '%6.2f' % (self.bfact)
        return linestr

    def _print_specific(self):
        """ PDBQT characters [68:79] """
        linestr =  ' ' * 2                      # [66:68]
        linestr += '%8.3f' % (self.charge)      # [68:76]
        linestr += ' ' * 1                      # [76:77]
        linestr += '%-2s' % (self.atype)       # [77:79]
        #linestr += '\n'
        return linestr






################################################################################
########### routine to find rings and return a dict of root member idx, CM, etc#
################################################################################
def ring_finder(filename):
	"""Create a dictionary objects with ring information"""
	#Setting a new object to be manipulated.
	mol = openbabel.OBMol()
	obConversion = openbabel.OBConversion()
	# imput and output be the same.
	obConversion.SetInAndOutFormats("sdf", "pdbqt")
	obConversion.ReadString(mol, filename)
	# Set the output variable
	out = []
	#iterate over the rings
	contador =0
	for r in mol.GetSSSR():
		if r.IsAromatic():
			contador+=1
			Vx = np.array([])
			Vy = np.array([])
			Vz = np.array([])
			mass = np.array([])
			for obatom in openbabel.OBMolAtomIter(mol):
				if r.IsInRing(obatom.GetIdx()):
					Vxcurrent=np.array([obatom.GetX()])
					Vx=np.concatenate((Vx,Vxcurrent))
					Vycurrent=np.array([obatom.GetY()])
                                       	Vy=np.concatenate((Vy,Vycurrent))
					Vzcurrent=np.array([obatom.GetZ()])
                                       	Vz=np.concatenate((Vz,Vzcurrent))
					mass_current=np.array([obatom.GetAtomicMass()])	
					mass=np.concatenate((mass,mass_current))
					root_serial=obatom.GetIdx()
			#return Aromatic center to the current ring.
			center = center_of_mass(Vx,Vy,Vz,mass)
			# Save information in a dictionary.
			element = {'ring_number' : contador,
			'root_atom': int(root_serial),
        		'Type': str(r.GetType()),
               		'center': center,
			'NumAtoms': mol.NumAtoms()}
			out.append(element)

	return out

def mod_pdbqt(inputfile,name,ring_list):
	"""Creates a list of PDBQT atom objects"""
	#abro el archivo a modificar
	f = inputfile.split("\n")
	# Defino la lista de anillos.
	#ring_list = ring_finder(inputfile)
	# si no encuentro aromaticos entonces devuelvo el mismo archivo.
	if len(ring_list)==0:
		# out variable
	        out = []
		print "WARNING: the input file ", name, " does not contain aromatic rings"
		for line in f:
			out.append(str(line) + "\n")
	else:
		current_dummy_serial = ring_list[0]['NumAtoms']
		#itero sobre los anillos.
		for element in ring_list:
			# out variable
	                out = []
			#itero sobre las lineas del archivo.
			for line in f:
				if line.startswith('ATOM  ') or line.startswith('HETATM'):
					#Defino un objeto parseable con las propiedades del atomo.
					atom = PDBQT(line)
					# defino el index del root del anillo.
					#the index for the root atom. O for furan, S for thiazole, N for pyrrole. For 6 membered aromatic rings, the first non carbon atom is used for root. For 5 members rings the O, S or N (BOSum=3, valence=3) will be used for root
					root_serial = element['root_atom']
					# Si encuentro el root en el archivo
					if int(root_serial)==int(atom.serial):
						# Defino el centro de masa el anillo.
						ring_CM  = element['center']
						#Genero una nueva linea para el CM.
						new_line = atom.getline()
						add_line = PDBQT(new_line) 
						# Modificar coordenadas con las del CM.
						current_dummy_serial = int(current_dummy_serial) + 1
						add_line.serial = current_dummy_serial
						add_line.x = ring_CM[0]
						add_line.y = ring_CM[1]
						add_line.z = ring_CM[2]
						# Modifico propiedades del CM.
						add_line.name = str('CM ')
						add_line.resNum = atom.resNum
						add_line.atype = 'E'
						add_line.charge = float(0.00)
						# Imprimo la linea root y CM.
						out.append(str(atom.getline()))
		                        	out.append(str(add_line.getline()))
					#imprimo las lineas hayan o no sido modificadas.
					else:				
						out.append(str(atom.getline()))
				else:
				#imprimo las lineas relacionadas con los REMARKS o torsiones ENDS etc.
					out.append(str(line))
			# defino el archivo f como el nuevo con las modificaciones.
			f=out
	return f

def add_dummy(input_filename, outfilename):
	# Defino el nombre del archivo y la extension.
	filename, file_extension = os.path.splitext(input_filename)
	if (file_extension==".sdf" or file_extension==".sd"):
		#molecule number
		N=0
		for mymol in pybel.readfile(file_extension.strip("."), input_filename):
			# Send file to a variable.
        		mol=mymol.write("pdbqt")
			# nombre de la molecula N.
			N+=1
			outfilename_base = mymol.title + "_" + str(N)
			#abro el archivo de salida
        		out = open(outfilename_base + '_dum.pdbqt', 'w')
        		# Hago las cosas.
			ring_list = ring_finder(mymol.write("sdf"))	
        		coordinates= mod_pdbqt(mol,outfilename_base + '_dum.pdbqt',ring_list)
        		# Escribo en disco.
        		out.writelines(["%s\n" % item  for item in coordinates[:-1]])
        		#out.write(coordinates)
        		out.close()
	if (file_extension==".pdb"):
		#defino el nombre de salida
        	salida = outfilename
        	#en caso de no haber dado archivo de salida seteamos el default
        	if salida is None:
                	salida = os.path.splitext(basename(filename))[0]+"_dum.pdbqt"
		for mymol in pybel.readfile(file_extension.strip("."), input_filename):
                        # Send file to a variable.
                        mol=mymol.write("pdbqt")
			print(mol)
                        # nombre de la molecula N.
                        #abro el archivo de salida
                        out = open(salida, 'w')
                        # Hago las cosas.
			ring_list = ring_finder(mymol.write("sdf"))
                        coordinates= mod_pdbqt(mol,salida,ring_list)
                        # Escribo en disco.
                        out.writelines(["%s\n" % item  for item in coordinates[:-1]])
                        #out.write(coordinates)
                        out.close()
	if file_extension==".pdbqt":
		#defino el nombre de salida
                salida = outfilename
                #en caso de no haber dado archivo de salida seteamos el default
                if salida is None:
                        salida = os.path.splitext(basename(filename))[0]+"_dum.pdbqt"
		with open(input_filename, 'r') as f:
    			fileContents = f.read() #read entire file
		mol=fileContents
		out = open(salida, 'w')
                # Hago las cosas.
		mymol = pybel.readfile("pdbqt",input_filename).next()
		ring_list = ring_finder(mymol.write("sdf"))
                coordinates= mod_pdbqt(mol,salida,ring_list)
                # Escribo en disco.
                out.writelines(["%s\n" % item  for item in coordinates[:-1]])
                #out.write(coordinates)
                out.close()



##############
#####MAIN##### 
##############



if __name__ == '__main__':
    import sys
    import getopt


    def usage():
        "Print helpful, accurate usage statement to stdout."
        print "Usage: prepare_dummies.py -l filename"
        print
        print "    Description of command..."
        print "         -l     ligand_filename (.pdbqt, .pdb or .sdf format)"
        print "    Optional parameters:"
        print "        [-o pdbqt_filename] (default output filename is ligand_filename_stem + _dum + .pdbqt)"
	print "        				In sdf files the output name will be taken from the title section.\n"
        print "							Developed by Elias Lopez (2017) contact:eliaslopez@qb.fcen.uba.ar\n"


    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'l:o:')
    except getopt.GetoptError, msg:
        print 'prepare_dummies.py: %s' %msg
        usage()
        sys.exit(2)

    # initialize required parameters
    #-l: ligand
    ligand_filename =  None
    # optional parameters
    #-o outputfilename
    outputfilename = None

    #'l:vo:d:A:CKU:B:R:MFI:Zgs'
    for o, a in opt_list:
        #print "o=", o, " a=", a
        if o in ('-l', '--l'):
            ligand_filename = a
        if o in ('-o', '--o'):
            outputfilename = a
        if o in ('-h', '--'):
            usage()
            sys.exit()
    if not  ligand_filename:
        print 'prepare_aromatic_center: ligand filename must be specified.'
        usage()
        sys.exit()
    else:
		# ejecuto el programa propiamente dicho.
		add_dummy(ligand_filename, outputfilename) 
		

