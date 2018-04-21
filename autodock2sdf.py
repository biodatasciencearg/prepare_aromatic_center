#!/usr/bin/python
import pybel as pybel
import openbabel
class Model:
	def __init__(self, model, clusterrank, clustersize, rmsdreference, freeenergybinding, ki, interenergy, vmdwhbonddesolv, totalinternalenergy, torsinalfreeenergy, unboundenergy, coordinates):
        	self.model = model
		self.ClusterRank = clusterrank
		self.ClusterSize = clustersize
		self.RmsdRefence = rmsdreference
		self.FreeEnergyBinding = freeenergybinding 
		self.Ki = ki
		self.InterEnergy = interenergy
		self.VdwHbondDesolv = vmdwhbonddesolv
		self.TotalInternalEnergy = totalinternalenergy
		self.TorsionalFreeEnergy = torsinalfreeenergy
		self.UnboundEnergy = unboundenergy
		self.coordinates = coordinates
	def SdfInfo(self):
        	linestr = ''
        	linestr += '\n>  <AUTOSCORE.MODEL>\n%-6i\n\n' % int(self.model)
		linestr += '>  <AUTOSCORE.RANK>\n%-6i\n\n'  % int(self.ClusterRank)
		linestr += '>  <AUTOSCORE.CLUSTERSIZE>\n%-6i\n\n'  % int(self.ClusterSize)
		linestr += '>  <AUTOSCORE.RMSD>\n%-10.3f\n\n' % float(self.RmsdRefence)
		linestr += '>  <AUTOSCORE.FreeEnergy>\n%-10.3f\n\n' % float(self.FreeEnergyBinding)
		linestr += '>  <AUTOSCORE.ki>\n%-10.3f\n\n' % float(self.Ki)
		linestr += '>  <AUTOSCORE.InterEnergy>\n%-10.3f\n\n' % float(self.InterEnergy)
		linestr += '>  <AUTOSCORE.VdwHbondDesolv>\n%-10.3f\n\n' % float(self.VdwHbondDesolv)
		linestr += '>  <AUTOSCORE.TotalInternalEnergy>\n%-10.3f\n\n' % float(self.TotalInternalEnergy)
		linestr += '>  <AUTOSCORE.TorsionalFreeEnergy>\n%-10.3f\n\n' % float(self.TorsionalFreeEnergy)
		linestr += '>  <AUTOSCORE.UnboundEnergy>\n%-10.3f\n\n' % float(self.UnboundEnergy)
		linestr += '$$$$\n'
		return linestr



def dlg_parser (pdbin):
	my_models = []
	s = ""
	with open(pdbin) as inf:
		canPrintLines = False # Boolean variable
		for line in inf:
    			if line.startswith("Keeping original residue number (specified in the input PDBQ file) for outputting."):
				canPrintLines = True
    			if canPrintLines:
				column = line.split() # split line into columns
				if line.startswith("MODEL"):
					model =  column[1]
				elif line.startswith("USER    Cluster Rank"):
					clusterrank= column[4]
				elif line.startswith("USER    Number of conformations in this cluster"):
					cluster_size = column[8]
				elif line.startswith("USER    RMSD from reference structure"):
					rmsd_reference = column[6]
				elif line.startswith("USER    Estimated Free Energy of Binding"):
					free_energy_binding = column[7]
				elif line.startswith("USER    Estimated Inhibition Constant, Ki"):
					Ki = column[6]
				elif line.startswith("USER    (1) Final Intermolecular Energy"):
					Inter_energy = column[6]
				elif line.startswith("USER        vdW + Hbond + desolv Energy"):
					VdwHbondDesolv = column[8]
				elif line.startswith("USER    (2) Final Total Internal Energy"):
					TotalInternalEnergy = column[7]
				elif line.startswith("USER    (3) Torsional Free Energy"):
					TorsionalFreeEnergy = column[6]
				elif line.startswith("USER    (4) Unbound System's Energy"):
					UnboundEnergy = column[7]
				elif line.startswith("ATOM"):
					s += "".join(list(line)[:70]) + "       " + "".join(list(line)[13]) + "\n"
				elif line.startswith("ENDMDL"):
					my_models.append(Model(model, clusterrank, cluster_size, rmsd_reference, free_energy_binding, Ki, Inter_energy, VdwHbondDesolv, TotalInternalEnergy, TorsionalFreeEnergy, UnboundEnergy, s))
					s = ""
	return my_models	

def dlg2sdf (dlgfile, name):
	s=""
	#itero sobre los modelos del pdb.
	for m in dlg_parser(dlgfile):
		mol=pybel.readstring("pdb", m.coordinates)
		mol.title = name
		mysdf=str(mol.write("sdf"))
		s += mysdf.strip("\n$$$$") + m.SdfInfo()
	return s
print(dlg2sdf("dock_out.pdb", "eliassoCra" ))


