# -*- coding: utf-8 -*-
import os, pickle, time, sys
from openmm import *
from openmm.app import *
from openmm.unit import kelvin, picosecond, femtosecond, bar, kilojoule_per_mole

WORK = r"C:\dockwork\md"
DT = 4.0  # fs
PROD_NS = float(sys.argv[1]) if len(sys.argv) > 1 else 10.0

with open(os.path.join(WORK, "modeller.pkl"), "rb") as f:
    blob = pickle.load(f)
topology, positions, n_prot = blob["topology"], blob["positions"], blob["n_prot"]
with open(os.path.join(WORK, "system.xml")) as f:
    system = XmlSerializer.deserialize(f.read())

# 生产用 NPT：先在体系里加 barostat，避免之后 reinitialize
system.addForce(MonteCarloBarostat(1.0 * bar, 300 * kelvin, 25))

# 蛋白重原子约束
heavy = [a.index for a in topology.atoms()
         if a.index < n_prot and a.element != Element.getBySymbol("H")]
rf = CustomExternalForce("0.5*k*periodicdistance(x,y,z,x0,y0,z0)^2")
rf.addGlobalParameter("k", 1000.0)
rf.addPerParticleParameter("x0"); rf.addPerParticleParameter("y0"); rf.addPerParticleParameter("z0")
for i in heavy:
    x, y, z = positions[i]
    rf.addParticle(i, [x, y, z])
system.addForce(rf)

integ = LangevinMiddleIntegrator(300 * kelvin, 1 / picosecond, DT * femtosecond)
sim = Simulation(topology, system, integ, Platform.getPlatformByName("OpenCL"))
sim.context.setPositions(positions)

dcd = os.path.join(WORK, "traj.dcd")
chk = os.path.join(WORK, "production.chk")
# 断点续跑
resume = os.path.exists(chk)
if resume:
    sim.loadCheckpoint(chk)
    done = sim.currentStep
    print(f"resumed from step {done}", flush=True)
    reporter = DCDReporter(dcd, 10000, append=True)
else:
    print("minimizing...", flush=True)
    LocalEnergyMinimizer.minimize(sim.context, 1000)
    sim.context.setVelocitiesToTemperature(300 * kelvin)
    done = 0
    reporter = DCDReporter(dcd, 10000)
sim.reporters.append(reporter)

# ---------- 平衡（仅新跑时）：NVT 50ps + NPT 200ps 逐步松约束 ----------
eq_steps = 0
if not resume:
    print("NVT equilibration 50 ps", flush=True)
    sim.step(12500)
    plan = [(500.0, 12500), (200.0, 12500), (50.0, 12500), (0.0, 12500)]
    for k, n in plan:
        sim.context.setParameter("k", k)
        sim.step(n)
    eq_steps = sim.currentStep
    print(f"equilibration done at step {eq_steps}", flush=True)
sim.context.setParameter("k", 0.0)

# ---------- 生产 ----------
target = eq_steps + int(PROD_NS * 1_000_000 / DT)
print(f"production target step {target} (~{PROD_NS} ns)", flush=True)
block = 100000
while sim.currentStep < target:
    t0 = time.time()
    sim.step(min(block, target - sim.currentStep))
    sim.saveCheckpoint(chk)
    pe = sim.context.getState(getEnergy=True).getPotentialEnergy()
    elapsed = time.time() - t0
    ns_done = sim.currentStep * DT / 1e6
    print(f"step {sim.currentStep} | {ns_done:.2f} ns | {pe} | {elapsed:.0f}s/block", flush=True)
print("PRODUCTION COMPLETE", flush=True)
