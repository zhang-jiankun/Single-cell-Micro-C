#!/usr/bin/env python 
#SBATCH -c 4
#SBATCH --mem 128G
#SBATCH --gpus 1 
#SBATCH --time 1-00:00:00
#SBATCH -J polychrom

import sys
import os
import numpy as np
import numpy.matlib
import h5py
import ast
import pandas as pd
import math

sys.path.append("/dshare/xielab/analysisdata/ThreeDGene/zhangjk/tools/polychrom-master")
sys.path.append("/dshare/xielab/analysisdata/ThreeDGene/zhangjk/tools/polychrom-master/Simulations")
from LEBondUpdater import bondUpdater

import polychrom
from polychrom.starting_conformations import grow_cubic
from polychrom.hdf5_format import HDF5Reporter, list_URIs, load_URI, load_hdf5_file
from polychrom.simulation import Simulation
from polychrom import polymerutils
from polychrom import forces
from polychrom import forcekits
import time

# ==== Load parameters
parValues = np.load(sys.argv[1], allow_pickle=True).item()

# general parameters 
trajectoryLength = parValues['trajectoryLength'] # time duration of simulation (down from 100,000)
density = parValues['density'] # density of the PBC box
  
#  ==========Extrusion sim parameters====================
# there is probably a more elegant way to read in text values than ast.literal_eval, but this works.  
GPU_ID = str(0)
N1 = parValues['monomers'] # Number of monomers in the polymer
M = parValues['replicates']  # concatinate replicates of polymer end-to-end (fewer TAD borders to assign, can be averaged if desired for LE)
N = N1 * M # number of monomers in the full simulation 
LIFETIME   = parValues['LEF_lifetime'] # 200 [Imakaev/Mirny use 200 as demo] extruder lifetime
SEPARATION = parValues['LEF_separation'] #  500 ave. separation between extruders in monomer units (extruder density) 

ctcfSites = parValues['ctcfSites'] # CTCF site locations [80,280,420,550,700,850] # positioned on HoxA
ctcfDir   = parValues['ctcfDir']
ctcfCapture = parValues['ctcfCapture'] # 0.9 80% capture probability per block if capture < than this, capture  
ctcfRelease = parValues['ctcfRelease'] # 0.003 % release probability per block. if capture < than this, release

interactionMatrix = parValues['interactionMatrix']
saveFolder = parValues['saveFolder']
saveFolder = saveFolder + '_NoLE' # important!

oneChainMonomerTypes = np.zeros(N1).astype(int) 

for site in parValues['EPSites']:
	oneChainMonomerTypes[site-3:site+3] = 1

# need to replicate and renormalize
loadProb = np.zeros(N1).astype(int) # discrete probability distribution that cohesin loads at site N
loadProb = loadProb + 0.0001
loadProb[parValues['EPSites'][-1]] = parValues['loadProb']
loadProb = numpy.matlib.repmat(loadProb, 1, M)
loadProb = loadProb/np.sum(loadProb) 

if not os.path.exists(saveFolder):
    os.mkdir(saveFolder)

lefPosFile = os.path.join(saveFolder, "LEFPos.h5")

# LEFNum = max(0,math.floor(N // SEPARATION )-1)
LEFNum = 0 # without loop extrusion

# less common parameters
attraction_radius = 1.5
num_chains = M  # simulation uses some equivalent chains  (5 in a real sim)
MDstepsPerCohesinStep = 800
smcBondWiggleDist = 0.2
smcBondDist = 0.5
saveEveryBlocks = 100   # save every 10 blocks (saving every block is now too much almost)
restartSimulationEveryBlocks = 100

# check that these loaded alright
print(f'LEF count: {LEFNum}')
print('interaction matrix:')
print(interactionMatrix)
print('monomer types:')
print(oneChainMonomerTypes)
print(saveFolder)
print('Starting simulation')


#==================================#
# Run 
#=================================#
import polychrom.lib.extrusion1Dv2 as ex1D # 1D classes 
ctcfLeftRelease = {}
ctcfRightRelease = {}
ctcfLeftCapture = {}
ctcfRightCapture = {}

# should modify this to allow directionality
for i in range(M): # loop over chains (this variable needs a better name Max)
    for t in range(len(ctcfSites)):
        pos = i * N1 + ctcfSites[t] 
        if ctcfDir[t] == 0:
            ctcfLeftCapture[pos] = ctcfCapture[t]  # if random [0,1] is less than this, capture
            ctcfLeftRelease[pos] = ctcfRelease[t]  # if random [0,1] is less than this, release
            ctcfRightCapture[pos] = ctcfCapture[t]
            ctcfRightRelease[pos] = ctcfRelease[t]
        elif ctcfDir[t] == 1: # stop Cohesin moving toward the right  
            ctcfLeftCapture[pos] = 0  
            ctcfLeftRelease[pos] = 1  
            ctcfRightCapture[pos] = ctcfCapture[t]
            ctcfRightRelease[pos] = ctcfRelease[t]
        elif ctcfDir[t] == 2:
            ctcfLeftCapture[pos] = ctcfCapture[t]  # if random [0,1] is less than this, capture
            ctcfLeftRelease[pos] = ctcfRelease[t]  # if random [0,1] is less than this, release
            ctcfRightCapture[pos] = 0
            ctcfRightRelease[pos] = 1
       
args = {}
args["ctcfRelease"] = {-1:ctcfLeftRelease, 1:ctcfRightRelease}
args["ctcfCapture"] = {-1:ctcfLeftCapture, 1:ctcfRightCapture}        
args["N"] = N 
args["LIFETIME"] = LIFETIME
args["LIFETIME_STALLED"] = LIFETIME  # no change in lifetime when stalled 

occupied = np.zeros(N)
occupied[0] = 1  # (I think this is just prevent the cohesin loading at the end by making it already occupied)
occupied[-1] = 1 # [-1] is "python" for end
cohesins = []

print('starting simulation with N LEFs=')
print(LEFNum)
for i in range(LEFNum):
    ex1D.loadOneFromDist(cohesins,occupied, args,loadProb) # load the cohesins 


with h5py.File(lefPosFile, mode='w') as myfile:
    
    dset = myfile.create_dataset("positions", 
                                 shape=(trajectoryLength, LEFNum, 2), 
                                 dtype=np.int32, 
                                 compression="gzip")
    steps = 100    # saving in 50 chunks because the whole trajectory may be large 
    bins = np.linspace(0, trajectoryLength, steps, dtype=int) # chunks boundaries 
    for st,end in zip(bins[:-1], bins[1:]):
        cur = []
        for i in range(st, end):
            ex1D.translocate(cohesins, occupied, args,loadProb)  # actual step of LEF dynamics 
            positions = [(cohesin.left.pos, cohesin.right.pos) for cohesin in cohesins]
            cur.append(positions)  # appending current positions to an array 
        cur = np.array(cur)  # when we finished a block of positions, save it to HDF5 
        dset[st:end] = cur
    myfile.attrs["N"] = N
    myfile.attrs["LEFNum"] = LEFNum
    
#=========== Load LEF simulation ===========#
trajectory_file = h5py.File(lefPosFile, mode='r')
LEFNum = trajectory_file.attrs["LEFNum"]  # number of LEFs
LEFpositions = trajectory_file["positions"]  # array of LEF positions  
steps = MDstepsPerCohesinStep # MD steps per step of cohesin  (set to ~800 in real sims)
Nframes = LEFpositions.shape[0] # length of the saved trajectory (>25000 in real sims)
print(f'Length of the saved trajectory: {Nframes}')
block = 0  # starting block 

# test some properties 
# assertions for easy managing code below 
assert (Nframes % restartSimulationEveryBlocks) == 0 
assert (restartSimulationEveryBlocks % saveEveryBlocks) == 0

savesPerSim = restartSimulationEveryBlocks // saveEveryBlocks
simInitsTotal  = (Nframes) // restartSimulationEveryBlocks
# concatinate monomers if needed
if len(oneChainMonomerTypes) != N:
    monomerTypes = np.tile(oneChainMonomerTypes, num_chains)
else:
    monomerTypes = oneChainMonomerTypes
    
N_chain = len(oneChainMonomerTypes)  
N = len(monomerTypes)
print(f'N_chain: {N_chain}')  # ~8000 in a real sim
print(f'N: {N}')   # ~40000 in a real sim
N_traj = trajectory_file.attrs["N"]
print(f'N_traj: {N_traj}')
assert N == trajectory_file.attrs["N"]
print(f'Nframes: {Nframes}')
print(f'simInitsTotal: {simInitsTotal}')
#==============================================================#
#                  RUN 3D simulation                              #
#==============================================================#
milker = bondUpdater(LEFpositions)
data = grow_cubic(N,int((N/(density*1.2))**0.333))  # starting conformation
reporter = HDF5Reporter(folder=saveFolder, max_data_length=50)
PBC_width = (N/density)**0.333
chains = [(N_chain*(k),N_chain*(k+1),0) for k in range(num_chains)]

for iteration in range(simInitsTotal):
    a = Simulation(N=N, 
                   error_tol=0.01, 
                   collision_rate=0.01, 
                   integrator ="variableLangevin", 
                   platform="cuda",
                   GPU = GPU_ID, 
                   PBCbox=(PBC_width, PBC_width, PBC_width),
                   reporters=[reporter],
                   precision="mixed")  # platform="CPU", 
    a.set_data(data)
    a.add_force(
        polychrom.forcekits.polymer_chains(
            a,
            chains=chains,
            nonbonded_force_func=polychrom.forces.heteropolymer_SSW,
            nonbonded_force_kwargs={
                'attractionEnergy': 0,  # base attraction energy for all monomers
                'attractionRadius': attraction_radius,
                'interactionMatrix': interactionMatrix,
                'monomerTypes': monomerTypes,
                'extraHardParticlesIdxs': []
            },
            bond_force_kwargs={
                'bondLength': 1,
                'bondWiggleDistance': 0.05
            },
            angle_force_kwargs={
                'k': 1.5
            }
        )
    )
    
    # ------------ initializing milker; adding bonds ---------
    # copied from addBond
    kbond = a.kbondScalingFactor / (smcBondWiggleDist ** 2)
    bondDist = smcBondDist * a.length_scale

    activeParams = {"length":bondDist,"k":kbond}
    inactiveParams = {"length":bondDist, "k":0}
    milker.setParams(activeParams, inactiveParams)
     
    # this step actually puts all bonds in and sets first bonds to be what they should be
    milker.setup(bondForce=a.force_dict['harmonic_bonds'],
                blocks=restartSimulationEveryBlocks)

    # If your simulation does not start, consider using energy minimization below
    if iteration == 0:
        a.local_energy_minimization() 
    else:
        a._apply_forces()
    
    for i in range(restartSimulationEveryBlocks):        
        if i % saveEveryBlocks == (saveEveryBlocks - 1):  
            a.do_block(steps=steps)
        else:
            a.integrator.step(steps)  # do steps without getting the positions from the GPU (faster)
        if i < restartSimulationEveryBlocks - 1: 
            curBonds, pastBonds = milker.step(a.context)  # this updates bonds. You can do something with bonds here
    data = a.get_data()  # save data and step, and delete the simulation
    del a
    
    reporter.blocks_only = True  # Write output hdf5-files only for blocks
    
    time.sleep(0.2)  # wait 200ms for sanity (to let garbage collector do its magic)

reporter.dump_data()
