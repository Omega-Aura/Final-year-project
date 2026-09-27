"""
Production MD with positional restraints on a named residue.

Written for the MAO-A system, where FAD is parameterised as its own GAFF2 residue rather than
as the covalently bound 8alpha-S-cysteinyl cofactor it actually is. Without the Cys linkage
nothing holds the flavin in place, so it is restrained to its crystallographic position. The
point of including FAD at all is that it forms one wall of the substrate cavity -- its absence
is what made every MAO-A docking before 2026-09-24 wrong -- and a restrained flavin reproduces
that wall even though the covalent bond is not modelled.

This is an approximation with a specific consequence: FAD cannot relax in response to the
ligand, so any induced fit involving the flavin is suppressed and the cofactor's own dynamics
are absent. For single-trajectory MM-GBSA of the LIGAND's binding energy that is defensible,
because FAD is part of the receptor on both sides of the subtraction and its internal energy
cancels. It is not a substitute for proper covalent parameterisation.

Protocol is otherwise identical to run_md_system.py so the MAO-A result stays comparable to
the TTBK1/TTBK2 runs.

Usage: python 06_md/run_md_restrained.py <system_dir> <restrained_resname> [steps] [k]
       k defaults to 10 kcal/mol/A^2 on heavy atoms only.
"""
import os
import sys
from openmm import app, unit, LangevinMiddleIntegrator, MonteCarloBarostat, Platform
from openmm import CustomExternalForce
from openmm.app import AmberPrmtopFile, AmberInpcrdFile, StateDataReporter, DCDReporter

SYS = sys.argv[1].rstrip("/\\")
RESNAME = sys.argv[2]
STEPS = int(sys.argv[3]) if len(sys.argv) > 3 else 5000000
K = float(sys.argv[4]) if len(sys.argv) > 4 else 10.0

p = lambda f: os.path.join(SYS, f)

prmtop = AmberPrmtopFile(p("complex.prmtop"))
inpcrd = AmberInpcrdFile(p("complex.inpcrd"))
system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0*unit.nanometer,
                             constraints=app.HBonds, rigidWater=True)

# Flat harmonic tether to the starting (crystal) coordinates. periodicdistance keeps this
# correct across the periodic boundary; a plain (x-x0)^2 would explode if the residue wrapped.
restraint = CustomExternalForce(
    "k*periodicdistance(x, y, z, x0, y0, z0)^2")
restraint.addGlobalParameter("k", K*unit.kilocalories_per_mole/unit.angstroms**2)
for name in ("x0", "y0", "z0"):
    restraint.addPerParticleParameter(name)

n = 0
for atom in prmtop.topology.atoms():
    if atom.residue.name == RESNAME and atom.element is not None \
            and atom.element.symbol != "H":
        restraint.addParticle(atom.index, inpcrd.positions[atom.index].value_in_unit(
            unit.nanometers))
        n += 1
if n == 0:
    sys.exit(f"ERROR: no heavy atoms found for residue '{RESNAME}' -- nothing restrained, "
             f"which would silently produce an unrestrained run")
system.addForce(restraint)
print(f"restrained {n} heavy atoms of residue {RESNAME} at k={K} kcal/mol/A^2", flush=True)

integrator = LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 2*unit.femtoseconds)
system.addForce(MonteCarloBarostat(1*unit.atmosphere, 300*unit.kelvin, 25))

platform = Platform.getPlatformByName('CUDA')
simulation = app.Simulation(prmtop.topology, system, integrator, platform)
simulation.context.setPositions(inpcrd.positions)
if inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

print(f'System {SYS}: {prmtop.topology.getNumAtoms()} atoms, {STEPS} production steps',
      flush=True)
print('Minimizing...', flush=True)
simulation.minimizeEnergy(maxIterations=2000)

print('Heating / equilibration (100 ps)...', flush=True)
simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
simulation.reporters.append(StateDataReporter(sys.stdout, 5000, step=True, temperature=True,
                                              potentialEnergy=True, volume=True, speed=True))
simulation.step(50000)

print('Production run...', flush=True)
simulation.reporters.append(DCDReporter(p('production.dcd'), 5000))
simulation.reporters.append(StateDataReporter(p('production.log'), 5000, step=True,
                                              time=True, potentialEnergy=True,
                                              temperature=True, volume=True, speed=True))
simulation.step(STEPS)

simulation.saveState(p('final_state.xml'))
positions = simulation.context.getState(getPositions=True).getPositions()
with open(p('final.pdb'), 'w') as f:
    app.PDBFile.writeFile(simulation.topology, positions, f)
print('Done.', flush=True)
