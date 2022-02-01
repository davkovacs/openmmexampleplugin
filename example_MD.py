"""
Script to create and run a hybrid ACE/MM simulation of a molecule in water. 
This uses the ACE custom openmm plugin
"""
from exampleplugin import ExampleForce # this is the ace plugin
from os import stat
from openmm import app, openmm, unit
from openmm.app import *
from openmm import *
from sys import stdout
from time import time
from openmmml.mlpotential import MLPotentialImplFactory, MLPotential, MLPotentialImpl
from typing import Iterable, Optional ,Iterable

# Use openmm-ml to define a new ML potenial and factory so we can use the createmixedsystem method for free
class ACEPoentialImpFratory(MLPotentialImplFactory):
    """This is a factory which creates the ACE potentials"""

    def createImpl(self, name: str, **args) -> MLPotentialImpl:
        return ACEPotentialImpl(name, **args)


class ACEPotentialImpl(MLPotentialImpl):
    """
    Impliment the ACE potential and add the force to a system.
    """

    def __init__(self, name: str, parameter_file: str) -> None:
        self.name = name
        self.parameter_file = parameter_file

    def addForces(self, topology: app.Topology, system: openmm.System, atoms: Optional[Iterable[int]], forceGroup: int, **args):
    
        # get the list of atoms to apply the ace force to
        includedAtoms = list(topology.atoms())
        if atoms is not None:
            includedAtoms = [includedAtoms[i] for i in atoms]
        # get the atomic numbers for each atom from the topology
        elements = [atom.element.atomic_number for atom in includedAtoms]

        # create the ace force  by loading the parameter file
        ace_force = ExampleForce(self.parameter_file)
        ace_force.setAtomicNumbers(elements)
        ace_force.setAtomInds(atoms)
        if topology.getPeriodicBoxVectors() is not None:
            ace_force.setUsesPeriodicBoundaryConditions(True)
        ace_force.setForceGroup(forceGroup)
        system.addForce(ace_force)

# now register to new potentials
MLPotential.registerImplFactory("ace", ACEPoentialImpFratory())


# running the simulation
ACE_PARAMETER_FILE = "../ACE1_bace17.json"   # add the file path here to the model 
# 1. build the normal mm system 
ace_potential = MLPotential("ace", parameter_file=ACE_PARAMETER_FILE)

pdb_file = app.PDBFile("../bace_water.pdb")
# work out the indices of the ligand atoms, we know they have id UNL
ace_atoms = [atom.index for atom in pdb_file.topology.atoms() if atom.residue.name == "UNL"]

# load the ligand and water force fields
ff = app.ForceField("../lig_CAT-17d.xml", "tip3p.xml")

# # create the normal system
mm_system = ff.createSystem(topology=pdb_file.topology, nonbondedMethod=app.PME, nonbondedCutoff=1 * unit.nanometer)

# # create the mixed system
hybrid_system = ace_potential.createMixedSystem(topology=pdb_file.topology, system=mm_system, atoms=ace_atoms)

temperature = 300 * unit.kelvin
hybrid_system.addForce(openmm.MonteCarloBarostat(1 * unit.atmosphere, temperature)) 
# build the simulation
time_step = 1 * unit.femtosecond
friction = 1 / unit.picosecond
temperature = 300 * unit.kelvin
integrator = openmm.LangevinIntegrator(temperature, friction, time_step)
platform = openmm.Platform.getPlatformByName("CPU")

simulation = app.Simulation(topology=pdb_file.topology, system=hybrid_system, integrator=integrator, platform=platform)
# set positions 
simulation.context.setPositions(pdb_file.positions)
# compute and print energies and forces
state = simulation.context.getState(getEnergy=True)
energy = state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
print(f"Starting system energy {energy}")
print("Optimizing the Geometry ...")
simulation.minimizeEnergy()

# compute properties again
state = simulation.context.getState(getEnergy=True, getPositions=True, enforcePeriodicBox=True)
energy = state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
print(f"minimised system energy {energy}")

# write out the minimised system to file
#Writer = app.dcdfile.DCDFile(file=open("BACE_NPT_MD.dcd"), topology=pdb_file.topology, dt=1*unit.femtosecond, interval=20)
#Writer.writeModel(positions=state.getPositions(), file=open("minimised_system.pdb", "w"))

data_freq=20
print("Running MD ...")
simulation.context.setVelocitiesToTemperature(250)
simulation.reporters.append(app.DCDReporter(f'BACE_NPT_MD.dcd', data_freq))
simulation.reporters.append(StateDataReporter(stdout, 50, step=True,
        potentialEnergy=True, temperature=True))
t0 = time() 
simulation.step(500000)
t1 = time()

print("500,000 step MD took {0:.2f} seconds".format(t1 - t0))
