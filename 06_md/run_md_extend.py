"""
Continue production MD from the saved state of the first 1 ns run, to tighten
the MM-GBSA SEM before reporting a result (see LOGBOOK.md 2026-08-21 entry).
"""
import sys
from openmm import app, unit, LangevinMiddleIntegrator, MonteCarloBarostat, Platform
from openmm.app import AmberPrmtopFile, DCDReporter, StateDataReporter

prmtop = AmberPrmtopFile('06_md/system/complex.prmtop')

system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0*unit.nanometer,
                              constraints=app.HBonds, rigidWater=True)

integrator = LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 2*unit.femtoseconds)
system.addForce(MonteCarloBarostat(1*unit.atmosphere, 300*unit.kelvin, 25))

platform = Platform.getPlatformByName('CUDA')
simulation = app.Simulation(prmtop.topology, system, integrator, platform)
simulation.loadState('06_md/system/final_state.xml')

simulation.reporters.append(DCDReporter('06_md/system/production2.dcd', 5000))
simulation.reporters.append(StateDataReporter(sys.stdout, 5000, step=True, time=True,
                                               potentialEnergy=True, temperature=True,
                                               volume=True, speed=True))
simulation.reporters.append(StateDataReporter('06_md/system/production2.log', 5000, step=True,
                                               time=True, potentialEnergy=True, temperature=True,
                                               volume=True, speed=True))

STEPS = int(sys.argv[1]) if len(sys.argv) > 1 else 4500000  # default +9 ns @ 2 fs
simulation.step(STEPS)

simulation.saveState('06_md/system/final_state2.xml')
print('Extension done.', flush=True)
