import rebound
import reboundx
import numpy as np
from rebound import InterruptiblePool

def getstr(dt):
    return "{:.0e}".format(dt)
def sim(par):
    run, mp, tau_e, e0, tmax, Nout, damping, integrator = par
    np.random.seed(run)
    theta = 2*np.pi*np.random.rand()
    sim = rebound.Simulation()
    sim.integrator = integrator
    sim.ri_ias15.epsilon = 0.
    sim.dt = 1.e-2
    sim.add(m=1.)
    
    sim.add(m=mp,a=1.,e=e0, theta=theta)
    sim.add(m=mp,a=2.,e=0.2)
    sim.move_to_com() # Moves to the center of momentum frame  
    
    if damping is True:
        rebx = reboundx.Extras(sim)
        rebx.add_modify_orbits_direct()
    
        tau_es = rebx.modify_orbits_direct.tau_e
        tau_es[1] = tau_e    
    
    times = np.logspace(0,np.log10(tmax),Nout)
    e1, e2, a1, a2 = np.zeros(len(times)), np.zeros(len(times)), np.zeros(len(times)), np.zeros(len(times))
    for i,time in enumerate(times):
        sim.integrate(time, exact_finish_time=1)
        o = sim.calculate_orbits()
        a1[i] = o[0].a
        a2[i] = o[1].a
        e1[i] = o[0].e
        e2[i] = o[1].e
        
    sim = rebound.Simulation()
    sim.integrator = "ias15"
    sim.ri_ias15.epsilon = 0.
    sim.dt = 1.e-3
    sim.add(m=1.)
    
    sim.add(m=mp,a=1.,e=e0, theta=theta)
    sim.add(m=mp,a=2.,e=0.2)
    sim.move_to_com() # Moves to the center of momentum frame  
    
    if damping is True:
        rebx = reboundx.Extras(sim)
        rebx.add_modify_orbits_direct()
    
        tau_es = rebx.modify_orbits_direct.tau_e
        tau_es[1] = tau_e    
    
    e1ref, e2ref, a1ref, a2ref = np.zeros(len(times)), np.zeros(len(times)), np.zeros(len(times)), np.zeros(len(times))
    for i,time in enumerate(times):
        sim.integrate(time, exact_finish_time=1)
        o = sim.calculate_orbits()
        a1ref[i] = o[0].a
        a2ref[i] = o[1].a
        e1ref[i] = o[0].e
        e2ref[i] = o[1].e

    a1err = np.fabs(a1-a1ref)/a1ref
    a2err = np.fabs(a2-a2ref)/a2ref
    e1err = np.fabs(e1-e1ref)/e1ref
    e2err = np.fabs(e2-e2ref)/e2ref
    
    return [times, a1, a2, e1, e2, a1err, a2err, e1err, e2err]

mp = 1.e-3
tau_e = -100
e0 = 1.e-1
tmax=1e6
Nout=1000
damping=True
integrator="whfast"
Nruns=20

params = [(run, mp, tau_e, e0, tmax, Nout, damping, integrator) for run in range(Nruns)]
pool = InterruptiblePool()
res = np.array(pool.map(sim, params))
    
ts = res[:,0,:]
a1 = res[:,1,:]
a2 = res[:,2,:]
e1 = res[:,3,:]
e2 = res[:,4,:]
a1err = res[:,5,:]
a2err = res[:,6,:]
e1err = res[:,7,:]
e2err = res[:,8,:]

import pickle

pickle.dump((Nruns, ts, a1, a2, e1, e2, a1err, a2err, e1err, e2err), open("test.p", "wb"))
