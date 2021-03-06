#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 21:54:00 2021

@author: iara
"""

import meep as mp
import numpy as num
#from meep.materials import SiO2
#from meep.materials import Si3N4
import matplotlib.pyplot as plt

#BS com o material que iremos usar nos chips e o comprimento de onda tbm 0.8
#
d = 0.075
default_d = 0.2 # branch separition in GDS
p=8 #numero de vezes que a simulacao ocorrera
I =[]

res = 25        # pixels/μm
three_d = False # 3d calculation?
  


gdsII_file = 'triBS.gds'
CELL_LAYER = 0
PORT1_LAYER = 1
PORT2_LAYER = 2
PORT3_LAYER = 3
PORT4_LAYER = 4
PORT5_LAYER = 5
PORT6_LAYER = 6
SOURCE_LAYER = 7
UPPER_BRANCH_LAYER = 31
LOWER_BRANCH_LAYER = 32
MEIO_BRANCH_LAYER = 33

t_oxide = 1.0
t_Si = 0.22
t_air = 0.78

dpml = 1
cell_thickness = dpml+t_oxide+t_Si+t_air+dpml

oxide = mp.Medium(epsilon=2.0736)
silicon= mp.Medium(epsilon=4.0804)

lcen = 0.795
fcen = 1/lcen
df = 0.2*fcen

cell_zmax = 0.5*cell_thickness if three_d else 0
cell_zmin = -0.5*cell_thickness if three_d else 0
si_zmax = 0.5*t_Si if three_d else 10
si_zmin = -0.5*t_Si if three_d else -10


for i in range(p):
   
    # read cell size, volumes for source region and flux monitors,
    # and coupler geometry from GDSII file
    upper_branch = mp.get_GDSII_prisms(silicon, gdsII_file, UPPER_BRANCH_LAYER, si_zmin, si_zmax)
    lower_branch = mp.get_GDSII_prisms(silicon, gdsII_file, LOWER_BRANCH_LAYER, si_zmin, si_zmax)
    meio_branch = mp.get_GDSII_prisms(silicon, gdsII_file, MEIO_BRANCH_LAYER, si_zmin, si_zmax)
    
    cell = mp.GDSII_vol(gdsII_file, CELL_LAYER, cell_zmin, cell_zmax)
    
    p1 = mp.GDSII_vol(gdsII_file, PORT1_LAYER, si_zmin, si_zmax)
    p2 = mp.GDSII_vol(gdsII_file, PORT2_LAYER, si_zmin, si_zmax)
    p3 = mp.GDSII_vol(gdsII_file, PORT3_LAYER, si_zmin, si_zmax)
    p4 = mp.GDSII_vol(gdsII_file, PORT4_LAYER, si_zmin, si_zmax)
    p5 = mp.GDSII_vol(gdsII_file, PORT5_LAYER, si_zmin, si_zmax)
    p6 = mp.GDSII_vol(gdsII_file, PORT6_LAYER, si_zmin, si_zmax)
    
    src_vol = mp.GDSII_vol(gdsII_file, SOURCE_LAYER, si_zmin, si_zmax)

    # displace upper and lower branches of coupler (as well as source and flux regions)
    
    if d != default_d:
        delta_y=d-default_d
        delta = mp.Vector3(y=delta_y)
        p1.center += delta
        p3.center -= delta
        p4.center += delta
        p6.center -= delta
        src_vol.center += delta
        cell.size += 2*delta
        for np in range(len(lower_branch)):
            lower_branch[np].center -= delta
            for nv in range(len(lower_branch[np].vertices)):
                lower_branch[np].vertices[nv] -= delta
        for np in range(len(upper_branch)):
            upper_branch[np].center += delta
            for nv in range(len(upper_branch[np].vertices)):
                upper_branch[np].vertices[nv] += delta
       
    
    geometry = upper_branch+meio_branch+lower_branch
    
    if three_d:
        oxide_center = mp.Vector3(z=-0.5*t_oxide)
        oxide_size = mp.Vector3(cell.size.x,cell.size.y,t_oxide)
        oxide_layer = [mp.Block(material=oxide, center=oxide_center, size=oxide_size)]
        geometry = geometry+oxide_layer
    
    sources = [mp.EigenModeSource(src=mp.GaussianSource(fcen,fwidth=df),
                                  volume=src_vol,
                                  direction = mp.X,
                                  eig_band=1,
                                  eig_parity=mp.NO_PARITY if three_d else mp.EVEN_Y+mp.ODD_Z,
                                  eig_match_freq=True)]
    
    sim = mp.Simulation(resolution=res,
                        cell_size=cell.size,
                        boundary_layers=[mp.PML(dpml)],
                        sources=sources,
                        geometry=geometry)
    
    mode1 = sim.add_mode_monitor(fcen, 0, 1, mp.ModeRegion(volume=p1))
    mode2 = sim.add_mode_monitor(fcen, 0, 1, mp.ModeRegion(volume=p2))
    mode3 = sim.add_mode_monitor(fcen, 0, 1, mp.ModeRegion(volume=p3))
    mode4 = sim.add_mode_monitor(fcen, 0, 1, mp.ModeRegion(volume=p4))
    mode5 = sim.add_mode_monitor(fcen, 0, 1, mp.ModeRegion(volume=p5))
    mode6 = sim.add_mode_monitor(fcen, 0, 1, mp.ModeRegion(volume=p6))
    
    sim.run(until_after_sources=100)
    
    # S parameters
    p1_coeff = sim.get_eigenmode_coefficients(mode1, [1], eig_parity=mp.NO_PARITY if three_d else mp.EVEN_Y+mp.ODD_Z).alpha[0,0,0]
    p2_coeff = sim.get_eigenmode_coefficients(mode2, [1], eig_parity=mp.NO_PARITY if three_d else mp.EVEN_Y+mp.ODD_Z).alpha[0,0,0]
    p3_coeff = sim.get_eigenmode_coefficients(mode3, [1], eig_parity=mp.NO_PARITY if three_d else mp.EVEN_Y+mp.ODD_Z).alpha[0,0,0 ]
    p4_coeff = sim.get_eigenmode_coefficients(mode4, [1], eig_parity=mp.NO_PARITY if three_d else mp.EVEN_Y+mp.ODD_Z).alpha[0,0,0]
    p5_coeff = sim.get_eigenmode_coefficients(mode5, [1], eig_parity=mp.NO_PARITY if three_d else mp.EVEN_Y+mp.ODD_Z).alpha[0,0,0]
    p6_coeff = sim.get_eigenmode_coefficients(mode6, [1], eig_parity=mp.NO_PARITY if three_d else mp.EVEN_Y+mp.ODD_Z).alpha[0,0,0]
    
    # transmittance
    #onde está a fonte
    norma=p2_coeff
    p1_trans = abs(p1_coeff)**2/abs()**2
    p2_trans = abs(p2_coeff)**2/abs(norma)**2
    p3_trans = abs(p3_coeff)**2/abs(norma)**2
    p4_trans = abs(p4_coeff)**2/abs(norma)**2
    p5_trans = abs(p5_coeff)**2/abs(norma)**2
    p6_trans = abs(p6_coeff)**2/abs(norma)**2
    


       # branch separation
    S =[]
    S.append(complex(d))
    S.append(p1_trans)
    S.append(p2_trans)
    S.append(p3_trans)
    S.append(p4_trans)
    S.append(p5_trans)
    S.append(p6_trans)
   
    I.append(S)
    d=d+0.03
    sim.reset_meep()
    
    
    
    
print(I)


# Data for plotting
I=num.transpose(I)
t = I[0]
S14=I[4]
S15 = I[5]
S16 = I[6]

fig, ax = plt.subplots()
ax.plot(t, S14,'r-o',label='S14')
ax.plot(t, S15,'g-o',label='S15')
ax.plot(t, S16,'y-o',label='S16')

plt.legend(loc = "upper right")

plt.xlabel("distances (μm)")

ax.grid()

fig.savefig("3BS_distances.png")
plt.show()

