#!/usr/bin/env python3

import numpy as np
import fileinput
import sys, string, os
import ctypes
import matplotlib.pyplot as plt
from lammps import lammps

# get dir name
dirname = os.path.join(os.getcwd(), "Feb_07")

# ---------------------------------------------------------
# step 1: read input data
#----------------------------------------------------------

# parse command line
argv = sys.argv
if len(argv) != 3:
  print(("Syntax: ./lammps_cyclo.py Npoints delta"))
  sys.exit()

# define number of points, and low and high box sizes
npoint = int(sys.argv[1])
delta = int(sys.argv[2])

# initial constants
hbar=0.658212           # Planck's constant (eV.fs/rad)
mub=5.78838e-2 # in meV/T
kb=8.61733e-2  # in meV/K
alat = 3.96

p_exchange = [-0.01768761631, 0.0, -1.98739396]
p_magelec = 0.18e-3
p_dmi = 0.06e-3

# expected minimum at ~620 Ang (exp. is 640 Ang)
emin = 620

# ---------------------------------------------------------
# step 2: run lammps on cell and output energy
#----------------------------------------------------------

def lammps_cyclo(cyclo_file,m_coeff,j_coeff):
  
  lmp = lammps(name='serial', cmdargs=["-screen","none"])
  # lmp = lammps()
  lmp.command("read_restart    %s" % cyclo_file)

  j1,j2,j3 = j_coeff[0],j_coeff[1],j_coeff[2]
  lmp.command("pair_style    hybrid/overlay spin/exchange 6.0 spin/magelec 4.5")
  lmp.command(("pair_coeff  * * spin/exchange exchange 6.0 %g %g %g"
      %(j1,j2,j3)))
  lmp.command(("pair_coeff  * * spin/magelec magelec 4.5 %g 0.0 0.0 1.0"
      %(m_coeff)))

  lmp.command("fix          1 all nve/spin lattice no")
  lmp.command("compute      out_mag all spin")
  lmp.command("variable     magx equal c_out_mag[1]")
  lmp.command("variable     magy equal c_out_mag[2]")
  lmp.command("variable     magz equal c_out_mag[3]")
  lmp.command("variable     nmag equal c_out_mag[4]")
  lmp.command("variable     emag equal c_out_mag[5]")
  lmp.command("thermo_style custom v_magx v_magy v_magz v_nmag v_emag pe press etotal")
  lmp.command("thermo       1")
  lmp.command("variable     ptot equal press")
  lmp.command("run          0")

  N = lmp.get_natoms()
  eng = lmp.extract_variable("emag","all",0)
  press = lmp.extract_variable("ptot","all",0)
  magx = lmp.extract_variable("magx","all",0)
  magy = lmp.extract_variable("magy","all",0)
  magz = lmp.extract_variable("magz","all",0)
  nmag = lmp.extract_variable("nmag","all",0)

  return eng/N, magx, magy, magz, nmag 


# ---------------------------------------------------------
# step 2: exchange only
#----------------------------------------------------------

eex = np.empty(npoint)

# exchange coefficients
j_coeff = np.empty(3)
j_coeff[0] = p_exchange[0]
j_coeff[1] = p_exchange[1]
j_coeff[2] = p_exchange[2]

# magelec coefficients
m_coeff = 0.0

print("Compute exchange energies")
for i in range(1, npoint+1):
  
  s1 = '../initial_spirals/restart_'
  stot=s1+str(i)
  
  out = lammps_cyclo(stot,m_coeff,j_coeff)

  eex[i-1] = out[0]

# ---------------------------------------------------------
# step 2: magneto-electric only
#----------------------------------------------------------

eme = np.empty(npoint)

# exchange coefficients
j_coeff = np.empty(3)
j_coeff[0] = 0.0 
j_coeff[1] = 0.0
j_coeff[2] = 1.0

# magelec coefficients
m_coeff = p_magelec

print("Compute ME energies")
for i in range(1, npoint+1):
  
  s1 = '../initial_spirals/restart_'
  stot=s1+str(i)
  
  out = lammps_cyclo(stot,m_coeff,j_coeff)

  eme[i-1] = out[0]

# ---------------------------------------------------------
# step 3: exchange + magneto-electric
#----------------------------------------------------------

eeme = np.empty(npoint)

# exchange coefficients
j_coeff = np.empty(3)
j_coeff[0] = p_exchange[0]
j_coeff[1] = p_exchange[1]
j_coeff[2] = p_exchange[2]

# magelec coefficients
m_coeff = p_magelec

print("Compute Ex. + ME energies")
for i in range(1, npoint+1):
  
  s1 = '../initial_spirals/restart_'
  stot=s1+str(i)
  
  out = lammps_cyclo(stot,m_coeff,j_coeff)

  eeme[i-1] = out[0]


# ---------------------------------------------------------
# step 4: plot results
#----------------------------------------------------------

# create dx table
dx = np.empty(npoint)
for i in range(1, npoint+1):
  isize = i*delta
  dx[i-1] = isize*alat

fig = plt.figure()

ax1 = plt.subplot(211)
ax1.set_xlabel('x ($\AA$)')
ax1.set_ylabel('Energy (eV / atom)')
ax1.plot(dx, eex, 'g--', label='Ex.')
ax1.plot(dx, eme, 'r--', label='ME')
ax1.plot(dx, eeme, 'b-', label='Ex. + ME')
ax1.axvline(x=emin, color = 'r', linestyle = '--')
ax1.legend(loc="upper right")

ax2 = plt.subplot(212)
ax2.set_xlabel('x ($\AA$)')
ax2.set_ylabel('Energy (quad)')
ax2.plot(dx, eeme, 'b-', label='Ex. + ME')
ax2.axvline(x=emin, color = 'r', linestyle = '--')
ax2.legend(loc="upper right")
ymin=-1.4705e-2
ymax=-1.47e-2
ax2.set_ylim((ymin, ymax)) 

plt.show()
fig.savefig(os.path.join(os.getcwd(), "test_cycloidal.pdf"), bbox_inches="tight")
plt.close(fig)
