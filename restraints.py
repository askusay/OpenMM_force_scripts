"""
OpenMM custom force scripts written as functions
https://github.com/askusay/OpenMM_force_scripts

See accompanying note-book for exmaples
"""

from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *

def print_rest(rest_type,**vars):
    print('Applying %s force:' % (rest_type))
    for key,value in vars.items():
        print('%s: %s = %d*%s' % (key, *value))

    return

forcedict = dict()

def pos_rest(system, pdb, sele, force, force_var):
    """Restrain atoms to the original locations in a PDB file
	Selection must be provide as list with atomic indicies"""

    assert len(sele) != 0, 'Atom list is empty!'

    # replace spaces with _
    force_var.replace(' ','_')

    print_rest('positional restraint',force_variable=(force_var, force, 'kilocalories_per_mole/angstroms**2'))
    forcedict[force_var]= '%s*%s' % (force, 'kilocalories_per_mole/angstroms**2')

    posforce = CustomExternalForce('%s*periodicdistance(x, y, z, x0, y0, z0)^2' % (force_var))
    posforce.addGlobalParameter(force_var, float(force)*kilocalories_per_mole/angstroms**2)
    posforce.addPerParticleParameter("x0")
    posforce.addPerParticleParameter("y0")
    posforce.addPerParticleParameter("z0")
    for i, z in zip(pdb.topology.atoms(), pdb.getPositions()):
        if i.index in sele:
            posforce.addParticle(int(i.index), [z[0], z[1], z[2]])
    system.addForce(posforce)

	
    return system

def z_rest(system, pdb, sele, force_var='z_force', force=1):
    """Restrain atoms to the original locations in a PDB file (Z-direction alone)
    Selection must be provide as list with atomic indicies"""
    
    assert len(sele) != 0, 'Atom list is empty!'

    # replace spaces with _
    force_var.replace(' ','_')

    print_rest('Z-direction restraint',force_variable=(force_var, force, 'kilocalories_per_mole/angstroms**2'))
    forcedict[force_var]= '%s*%s' % (force, 'kilocalories_per_mole/angstroms**2')

    boxlz = system.getDefaultPeriodicBoxVectors()[2][2].value_in_unit(nanometer)
    
    z_res = CustomExternalForce('%s*(zabs^2); \
                                    zabs = min(zd, boxlz-zd); \
                                    zd = abs(z-z0)' % (force_var))
    z_res.addGlobalParameter('boxlz', boxlz)
    z_res.addGlobalParameter(force_var,float(force)*kilocalories_per_mole/angstroms**2)
    z_res.addPerParticleParameter("z0")
    for i, z in zip(pdb.topology.atoms(), pdb.getPositions()):
        if i.index in sele:
            z_res.addParticle(int(i.index), [z[2]])
            
    system.addForce(z_res)
    
    return system

def fb_atoms_rest(system, pdb, sele, dist_var, force_var, distance, force=5.92):
    """Flat-bottom restrains between 2 atoms
    Selection must be a list containing 2 atoms or lists of lists containing 2 atoms each:
    [a1,a2] or [[a1,a2],[a3,a4]...]
    Please note that selections containing multiple selections will be assigned the same dist_var and force_var
    """

    assert type(sele) == list, 'Selection must be provided in a list'
    if type(sele[0]) == list:
        assert all([len(i) == 2 for i in sele]), 'All selections must contain 2 atoms'
    else:
        assert len(sele) == 2, 'Selection must contain 2 atoms'
        #encase single selection in a list
        sele = [sele]

    # replace spaces with _
    dist_var.replace(' ','_')
    force_var.replace(' ','_')

    print_rest('Flat-bottom atom restraint',force_variable=(force_var, force, 'kilocalories_per_mole/angstroms**2'), distance_variable=(dist_var, distance, 'Angstrom'))
    forcedict[force_var]= '%s*%s' % (force, 'kilocalories_per_mole/angstroms**2')
    forcedict[dist_var]= '%s*%s' % (distance, 'Angstrom')

    fb_atoms = CustomBondForce('step(r-%s) * (%s/2) * (r-%s)^2' % (dist_var,    force_var,dist_var))
    fb_atoms.addGlobalParameter(dist_var,float(distance)*angstroms)
    fb_atoms.addGlobalParameter(force_var,float(force)*kilocalories_per_mole/angstroms**2)

    for s in sele:
        a1, a2 = s
        fb_atoms.addBond(a1,a2)
		
    fb_atoms.setUsesPeriodicBoundaryConditions(True)
    system.addForce(fb_atoms)

    return system
 
def fb_groups_rest(system, pdb, sele, dist_var, force_var, distance, force=5.92):
    """Flat-bottom restrains between 2 groups of atoms
    Selection must be a list containing 2 lists of atoms:
    [[a1,a2,a3],[a4,a5,a6]]"""

    assert all([type(s)==list for s in sele]), 'Selection must be a list containing 2 lists'

    # replace spaces with _
    force_var.replace(' ','_')
    force_var.replace(' ','_')

    print_rest('Flat-bottom centroid restraint',force_variable=(force_var, force, 'kilocalories_per_mole/angstroms**2'), distance_variable=(dist_var, distance, 'Angstrom'))
    forcedict[force_var]= '%s*%s' % (force, 'kilocalories_per_mole/angstroms**2')
    forcedict[dist_var]= '%s*%s' % (distance, 'Angstrom')
    
    fb_groups = CustomCentroidBondForce(2, 'step(distance(g1,g2)-%s) * (%s/2)*(distance(g1,g2)-%s)^2' % (dist_var,force_var,dist_var))
    fb_groups.addGlobalParameter(dist_var,float(distance)*angstroms)
    fb_groups.addGlobalParameter(force_var,float(force)*kilocalories_per_mole/angstroms**2)
    fb_groups.addGroup(sele[0])
    fb_groups.addGroup(sele[1])
    bondGroups = [0, 1]
    fb_groups.addBond(bondGroups,[])
    fb_groups.setUsesPeriodicBoundaryConditions(True)
    system.addForce(fb_groups)

    return system

def angle_rest(system,pdb,sele,angle_var,force_var,angle,force):
    """Angle restraint between 3 atoms
    Selection must be a list containing 3 atoms or lists of lists containing 3 atoms each:
    [a1,a2,a3] or [[a1,a2,a3],[a4,a5,a6]...]
    Please note that selections containing multiple selections will be assigned the same angle_var and force_var"""
	
    assert type(sele) == list, 'Selection must be provided in a list'
    if type(sele[0]) == list:
        assert all([len(i) == 3 for i in sele]), 'All selections must contain 3 atoms'
    else:
        assert len(sele) == 3, 'Selection must contain 3 atoms'

    # replace spaces with _
    force_var.replace(' ','_')
    angle_var.replace(' ','_')

    print_rest('Angle restraint',force_variable=(force_var, force, 'kilocalories_per_mole/radian**2'), angle_variable=(angle_var, angle, 'radians'))
    forcedict[force_var]= '%s*%s' % (force, 'kilocalories_per_mole/radian**2')
    forcedict[angle_var]= '%s*%s' % (angle, 'radians')

    angle_f = CustomAngleForce("0.5*%s*(theta-%s)^2" % (force_var,angle_var))
    angle_f.addGlobalParameter(angle_var, float(angle)*radians)
    angle_f.addGlobalParameter(force_var, float(force)*kilocalories_per_mole/angstroms**2)

    for s in sele:
        a1,a2,a3 = s
        angle_f.addAngle(a1,a2,a3)
	
    angle_f.setUsesPeriodicBoundaryConditions(True)
    system.addForce(angle_f)

    return system
	
def torsion_rest(system,pdb,sele,angle_var,force_var,angle,force):
    """Torsional restraint between 4 atoms
    Selection must be a list containing 3 atoms or lists of lists containing 3 atoms each:
    [a1,a2,a3,a4] or [[a1,a2,a3,a4],[a5,a6,a7,a8]...]
    Please note that selections containing multiple selections will be assigned the same angle_var and force_var"""

    assert type(sele) == list, 'Selection must be provided in a list'
    if type(sele[0]) == list:
        assert all([len(i) == 4 for i in sele]), 'All selections must contain 4 atoms'
    else:
        assert len(sele) == 4, 'Selection must contain 4 atoms'
        #encase single selection in a list
        sele = [sele]

    # replace spaces with _
    force_var.replace(' ','_')
    angle_var.replace(' ','_')

    print_rest('Torsional restraint',force_variable=(force_var, force, 'kilocalories_per_mole'), angle_variable=(angle_var, angle, 'radians'))
    forcedict[force_var]= '%s*%s' % (force, 'kilocalories_per_mole')
    forcedict[angle_var]= '%s*%s' % (angle, 'radians')

    torsion = CustomTorsionForce("0.5*%s*(1-cos(theta-%s))" % (force_var,angle_var))
    torsion.addGlobalParameter(angle_var, float(angle)*radians)
    torsion.addGlobalParameter(force_var, float(force)*kilocalories_per_mole)

    for s in sele:
        a1,a2,a3,a4 = s
        torsion.addTorsion(a1,a2,a3,a4)
		
    torsion.setUsesPeriodicBoundaryConditions(True)
    system.addForce(torsion)

    return system

def forces():
    return forcedict