#!/usr/bin/env python3

import sys, string, os
import numpy as np
from lammps import lammps

# get dir name
dirname = os.path.join(os.getcwd(), "Nov_13")

# ---------------------------------------------------------
# step 1: read input data
#----------------------------------------------------------

# parse command line
argv = sys.argv
if len(argv) != 3:
  print(("Syntax: ./create_spirals.py Npoints delta"))
  sys.exit()

# define number of points, and low and high box sizes
npoint = int(sys.argv[1])
delta = int(sys.argv[2])

# initial constants
mub=5.78838e-2 # in meV/T
kb=8.61733e-2  # in meV/K
alat = 3.96

# ---------------------------------------------------------
# step 2: define function generating cells
#----------------------------------------------------------

def lammps_spiral(i_size,index):
  
  # define supercell as 2x spiral size
  # this is a crude approx., to be improved
  Nx=int(i_size)
  Ny=2
  Nz=2
  
  # define number of spins in bcc cell
  N=Nx*Ny*Nz
  
  Lx=alat*Nx
  Ly=alat*Ny
  Lz=alat*Nz
  sx0=1.0
  sy0=0.0
  sz0=0.0
  
  # print(Lx,Ly,Lz)
  # print(N)

  lmp = lammps(name='serial', cmdargs=["-screen","none"])
  lmp.command("units 	    metal")
  lmp.command("atom_style   spin")
  lmp.command("dimension    3")
  lmp.command("boundary     p p p")
  lmp.command("atom_modify  map array")
  lmp.command("region 	    box block 0.0 %g 0.0 %g 0.0 %g" % (Lx,Ly,Lz))
  lmp.command("create_box   1 box")
 
  # create atoms
  p = 1
  for i in range(0, Nx):
    lx1 = alat*i 
    for j in range(0, Ny):
      ly1 = alat*j
      for k in range(0, Nz):
        lz1 = alat*k 
        x1 = lx1
        y1 = ly1
        z1 = lz1
        lmp.command("create_atoms 1 single %g %g %g" % (x1,y1,z1))
        p += 1
  
  # set spins
  p = 1
  for i in range(0, Nx):
    lx1 = alat*i 
    for j in range(0, Ny):
      ly1 = alat*j
      for k in range(0, Nz):
        lz1 = alat*k 
        x1 = lx1
        y1 = ly1
        z1 = lz1
        dot1 = 2.0*np.pi*x1/(Nx*alat)
        s1x=(sx0*np.cos(dot1))*(-1.0)**(i+j+k)
        s1y=0.0
        s1z=(sx0*np.sin(dot1))*(-1.0)**(i+j+k)
        isnorm1=1.0/(s1x*s1x+s1y*s1y+s1z*s1z)**0.5
        s1x*=isnorm1
        s1y*=isnorm1
        s1z*=isnorm1
        lmp.command("set atom %d spin 2.0 %g %g %g" % (p,s1x,s1y,s1z))
        p += 1

  # set masses
  lmp.command("mass     1 55.845")

  # output dump and restart files
  lmp.command("compute 	     outsp all property/atom spx spy spz sp fmx fmy fmz")
  lmp.command("dump          1 all custom 1 dump_%d type x y z c_outsp[1] c_outsp[2] c_outsp[3]" % index)

  lmp.command("run           0")
  lmp.command("write_restart restart_%d" % index)

# ---------------------------------------------------------
# step 3: loop generating lammps restart files
#----------------------------------------------------------

print("Generate cycloid inputs")
for i in range(1, npoint+1):
  
  isize = i*delta

  # writing lammps data file
    
  lammps_spiral(isize,i)

