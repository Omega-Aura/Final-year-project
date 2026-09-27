"""
Minimize, equilibrate and run production MD for any prepared complex.

Same protocol as run_md.py (ff14SB/gaff2/tip3p, PME 1.0 nm, HBonds constraints, Langevin 300 K
2 fs, Monte Carlo barostat 1 atm) -- the only change is that the system directory is an
argument instead of being hardcoded to 06_md/system, so TTBK2 and any later complex are run
under conditions identical to the TTBK1 baseline rather than a re-typed approximation of it.

Production is run in one continuous block. The TTBK1 baseline reached the same 10 ns as
1 ns + a 9 ns restart from a saved state, which is equivalent: same integrator, same barostat,
continuous dynamics.

Usage: python 06_md/run_md_system.py <system_dir> [production_steps]
       default 5000000 steps = 10 ns at 2 fs
"""
import os
import sys
from openmm import app, unit, LangevinMiddleIntegrator, MonteCarloBarostat, Platform
from openmm.app import AmberPrmtopFile, AmberInpcrdFile, StateDataReporter, DCDReporter

SYS = sys.argv[1].rstrip("/\\")
STEPS = int(sys.argv[2]) if len(sys.argv) > 2 else 5000000   # 10 ns @ 2 fs

p = lambda f: os.path.join(SYS, f)

prmtop = AmberPrmtopFile(p("complex.prmtop"))
inpcrd = AmberInpcrdFile(p("complex.inpcrd"))

system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0*unit.nanometer,
                             constraints=app.HBonds, rigidWater=True)

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
simulation.reporters.append(DCDReporter(p('production.dcd'), 5000))          # every 10 ps
simulation.reporters.append(StateDataReporter(p('production.log'), 5000, step=True,
                                              time=True, potentialEnergy=True,
                                              temperature=True, volume=True, speed=True))
simulation.step(STEPS)

simulation.saveState(p('final_state.xml'))
positions = simulation.context.getState(getPositions=True).getPositions()
with open(p('final.pdb'), 'w') as f:
    app.PDBFile.writeFile(simulation.topology, positions, f)
print('Done.', flush=True)
