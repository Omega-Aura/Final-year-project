"""
Minimize, equilibrate, and run production MD for cand_003 / TTBK1 (7JXX) complex.
Preliminary Week-2 run, time-boxed for the 2-day review window (see LOGBOOK.md
2026-08-21 entries) -- not the master plan's full multi-replicate scope.
"""
import sys
from openmm import app, unit, LangevinMiddleIntegrator, MonteCarloBarostat, Platform
from openmm.app import AmberPrmtopFile, AmberInpcrdFile, PDBReporter, StateDataReporter, DCDReporter

prmtop = AmberPrmtopFile('06_md/systems/system/complex.prmtop')
inpcrd = AmberInpcrdFile('06_md/systems/system/complex.inpcrd')

system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0*unit.nanometer,
                              constraints=app.HBonds, rigidWater=True)

integrator = LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 2*unit.femtoseconds)
system.addForce(MonteCarloBarostat(1*unit.atmosphere, 300*unit.kelvin, 25))

platform = Platform.getPlatformByName('CUDA')
simulation = app.Simulation(prmtop.topology, system, integrator, platform)
simulation.context.setPositions(inpcrd.positions)
if inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

print('Minimizing...', flush=True)
simulation.minimizeEnergy(maxIterations=2000)

print('Heating / NVT equilibration (100 ps)...', flush=True)
simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
simulation.reporters.append(StateDataReporter(sys.stdout, 5000, step=True, temperature=True,
                                               potentialEnergy=True, volume=True, speed=True))
simulation.step(50000)  # 100 ps NVT-ish equilibration (barostat active throughout, standard NPT equil practice)

print('Production run...', flush=True)
simulation.reporters.append(DCDReporter('06_md/systems/system/production.dcd', 5000))  # every 10 ps
simulation.reporters.append(StateDataReporter('06_md/systems/system/production.log', 5000, step=True,
                                               time=True, potentialEnergy=True, temperature=True,
                                               volume=True, speed=True))

PRODUCTION_STEPS = int(sys.argv[1]) if len(sys.argv) > 1 else 500000  # default 1 ns @ 2 fs
simulation.step(PRODUCTION_STEPS)

simulation.saveState('06_md/systems/system/final_state.xml')
positions = simulation.context.getState(getPositions=True).getPositions()
with open('06_md/systems/system/final.pdb', 'w') as f:
    app.PDBFile.writeFile(simulation.topology, positions, f)
print('Done.', flush=True)
