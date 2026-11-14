# -*- coding: utf-8 -*-
import os, pickle
from openmm import *
from openmm.app import *
from openmm.unit import nanometer, kelvin, picosecond, femtosecond, bar, kilojoule_per_mole

WORK = r"C:\dockwork\md"
with open(os.path.join(WORK, "modeller.pkl"), "rb") as f:
    blob = pickle.load(f)
topology, positions = blob["topology"], blob["positions"]
n_prot = blob["n_prot"]
with open(os.path.join(WORK, "system.xml")) as f:
    system = XmlSerializer.deserialize(f.read())

# 蛋白重原子索引（用于约束）
heavy = [a.index for a in topology.atoms()
         if a.index < n_prot and a.element != Element.getBySymbol("H")]
print("restrained protein heavy atoms:", len(heavy))

def make_restraint(k):
    rf = CustomExternalForce("0.5*k*periodicdistance(x,y,z,x0,y0,z0)^2")
    rf.addGlobalParameter("k", k)
    rf.addPerParticleParameter("x0"); rf.addPerParticleParameter("y0"); rf.addPerParticleParameter("z0")
    for i in heavy:
        x, y, z = positions[i]
        rf.addParticle(i, [x, y, z])
    return rf

# ---------- 最小化（强约束） ----------
rf = make_restraint(1000.0)
system.addForce(rf)
integ = LangevinMiddleIntegrator(300 * kelvin, 1 / picosecond, 4.0 * femtosecond)
sim = Simulation(topology, system, integ, Platform.getPlatformByName("OpenCL"))
sim.context.setPositions(positions)
print("minimizing...")
LocalEnergyMinimizer.minimize(sim.context, 1000)
e0 = sim.context.getState(getEnergy=True).getPotentialEnergy()
print("potential after minimization:", e0)

# ---------- NVT 50 ps（强约束） ----------
sim.context.setVelocitiesToTemperature(300 * kelvin)
sim.step(12500)   # 50 ps
print("NVT done, T:", sim.context.getState(getEnergy=True).getKineticEnergy())

# ---------- NPT 200 ps，逐步放松约束 ----------
system.addForce(MonteCarloBarostat(1.0 * bar, 300 * kelvin, 25))
# 重新建 context 以纳入 barostat
sim.context.reinitialize(preserveState=True)
for kval, nsteps in [(500.0, 12500), (200.0, 12500), (50.0, 12500), (0.0, 12500)]:
    sim.context.setParameter("k", kval)
    sim.step(nsteps)
st = sim.context.getState(getEnergy=True, getPositions=True, getVelocities=True)
print("NPT done, potential:", st.getPotentialEnergy())
print("box:", st.getPeriodicBoxVectors())

sim.saveCheckpoint(os.path.join(WORK, "equilibrated.chk"))
with open(os.path.join(WORK, "equilibrated.pkl"), "wb") as f:
    pickle.dump({"positions": st.getPositions(), "velocities": st.getVelocities(),
                 "box": st.getPeriodicBoxVectors()}, f)
print("saved equilibrated.chk")
