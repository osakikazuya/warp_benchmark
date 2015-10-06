# ioad warp and various scripts 
from warp           import *   # warp code 
from lattice        import *   # lattice definition scripts (MAD syntax)
from species        import *   # species handling for loading routines 
from errorcheck     import *   # contains scripts to check for input errors
from realboundaries import *   # implements conducting pipe boundary
                               #   conditions in field solve
from timedependentvoltage import * 
from gridcrossingdiags import *
from extpart import *
# load modules necessary to use PW and PR.
from PWpickle import *
from PRpickle import *

random.seed(100)               # set a seed of random numbers 
                               # to get the exactly same distribution.

# ---------- set derived and warp variables
np = 4 #particle number
spi = Species(type=Uranium, charge_state = +34, name = "Ion")
spi.vbeam=VELOSITY;
top.aion          =  spi.mass/top.amu
top.zion          =  spi.charge/top.echarge
top.emitx         =  0.                    # ** arb (particles accumulated) ** beam emittance (m-rad) 
top.npmax         =  np                 # simulation particle number 
top.lrelativ      =  true           # turn off relativity 
if (top.lrelativ):         # turn on relativistic self-field corrections
   top.relativity = 1
print spi.charge, spi.sw, spi.mass/top.emass


# ---------- define simulation box.
w3d.l2symtry = false      # 2-fold transverse symmetry
w3d.l4symtry = true      # 4-fold transverse symmetry
checksymmetry()           # check that inputs are consistent with symmetries
                          # using errorcheck package function.

# boundary condition
w3d.boundxy  = 0   
w3d.bound0   = 0
w3d.boundnz  = 0

# read data file of magnetic field
data_file = "B_table_R0.08_L0.3_B_DATA.dat"
data = getdatafromtextfile(data_file,nskip=0,dims=[4,None],)
[d_r,d_z,br,bz] = data[0:4]


bz00 = 1.7647058823529411 # Normalize by central magnet B_z on-axis.  Note: this is NOT the peak nonlinear field! 
br = br*BZFIELD/bz00
bz = bz*BZFIELD/bz00

# Calculate dr,dz for uniform mesh data   
d_dr = diff(unique(d_r))
dr = average(d_dr)

d_dz = diff(unique(d_z))
dz = average(d_dz)


# Generate data mesh and field arrays with +z and -z (reflected) coordinates 
# excluding boundary points.  df_ prefix for "data full".
#    Note:  z coordinates flipped in sign since d_z contains data from z = 0 to z = max (downstream) 
#           br data flipped in sign consistently 

f_z   = concatenate(( -d_z[-1:0:-1], d_z[:]))
f_r   = concatenate((  d_r[-1:0:-1], d_r[:]))
f_br  = concatenate((-br[-1:0:-1],br[:]))
f_bz  = concatenate(( bz[-1:0:-1],bz[:]))

#   r_m = radial mesh vector 
#   z_m = axial  mesh vector 

nr = nint((f_r.max()-f_r.min())/dr)
r_mmin = f_r.min()
r_mmax = f_r.max()

r_m = r_mmin + dr*arange(0,nr+1,dtype=float)

nz = nint((f_z.max()-f_z.min())/dz)
z_mmin = f_z.min()
z_mmax = f_z.max()

z_m = z_mmin+ dz*arange(0,nz+1,dtype=float)

# Load field data in arrays 
br_m = fzeros((nr+1,nz+1))
bz_m = fzeros((nr+1,nz+1))

ir_m = nint((f_r - r_mmin)/dr)
iz_m = nint((f_z - z_mmin)/dz)

ii = ir_m + iz_m*(nr+1)

br_m.ravel('F').put(ii,f_br)
bz_m.ravel('F').put(ii,f_bz)

# Save linear field data to an external binary array 
fo = PW("s4.rz.pkl")
fo.s4_dz = dz
fo.s4_dr = dr
fo.s4_nr = nr
fo.s4_nz = nz
fo.s4_r_mmin = r_mmin
fo.s4_r_mmax = r_mmax
fo.s4_z_mmin = z_mmin
fo.s4_z_mmax = z_mmax
fo.s4_r_m    = r_m
fo.s4_z_m    = z_m
fo.s4_br_m   = br_m
fo.s4_bz_m   = bz_m

# Add b_field to grids
so_br_m = fzeros((nr+1,1,nz+1))
so_br_m[:,0,:] = br_m
so_bz_m = fzeros((nr+1,1,nz+1))
so_bz_m[:,0,:] = bz_m

addnewbgrd(zs=min(z_m),ze=max(z_m),xs=0.,dx=dr,dy=1.,bx=so_br_m,bz=so_bz_m,
           rz=true,nx=fo.s4_nr,nz=fo.s4_nz)

# Grid bounds
lx_grid = 10*r_mmax      # length edge of simulation grid
lz_grid = 2*z_mmax
w3d.xmmin =  0.                # r-grid min limit
w3d.xmmax =  lx_grid/2.        # r-grid max limit
w3d.zmmin = -lz_grid/2*1.01    # z-grid min limit
w3d.zmmax =  lz_grid/2         # z-grid max limit

top.prwall = lx_grid/2

# Grid size for space charge calculation. It is not needed this time.
ngridxy = GRIDR
ngridz  = GRIDZ
symx = 1
symy = 1
if (w3d.l4symtry):
 symx = 0.5
 symy = 0.5
elif (w3d.l2symtry):
 symx = 0.5 

w3d.nx = int(symx*ngridxy) 
w3d.ny = int(symy*ngridxy)
w3d.nz = ngridz

nx1 = ngridxy
xmesh = w3d.xmmin + (w3d.xmmax-w3d.xmmin)*array(range(0,nx1+1))/nx1
if w3d.l4symtry:
 nx1 = ngridxy/2
 xmesh = w3d.xmmax*array(range(0,nx1+1))/nx1

ny1 = ngridxy
ymesh = w3d.ymmin + (w3d.ymmax-w3d.ymmin)*array(range(0,ny1+1))/ny1
if w3d.l4symtry or w3d.l2symtry:
 ny1 = ngridxy/2
 ymesh = w3d.ymmax*array(range(0,ny1+1))/ny1

nz1 = ngridz
zmesh = w3d.zmmin + (w3d.zmmax-w3d.zmmin)*array(range(0,nz1+1))/nz1

# create test particles coordinate array and some personal variables.
nstep    = NSTEP      # number s-steps per lattice period
nadvance = NADVANCE

# setting initial particle condition.
x  = array([4e-2] * np)
y  = array([0] * np)
vx  = array([0] * np)
vy  = array([0,3e2,3e3,3e4,])
z  = array([-1*z_mmax] * np)
vz  = array([VELOSITY] *np)

# output data to a pdb file.
ff=PW("out.pdb")
ff.np = np
ff.x = x
ff.y = y
ff.z = z
ff.vx = vx
ff.vy = vy
ff.vz = vz
ff.close()

# input data from a pdb file.
ff=PR("out.pdb")

# ---------- setup function to add particles
def addparts():
    if top.it==1 :
      print "Adding particles."
      #if ff.np != N_sim:
      #  raise "ff.np is not equal to N_sim."
      spi.addpart(x=ff.x,y=ff.y,z=ff.z,vx=ff.vx,vy=ff.vy,vz=ff.vz,gi=ones(ff.np),lallindomain=0)
      if top.it==1e6: 
        ff.close()
    return
installuserinjection(addparts)

#set simulation step size
top.dt            =  ((w3d.zmmax-w3d.zmmin)/float(nstep))/spi.vbeam
#calculation of magnetic pushing force
top.ibpush   = 2 # 0:not cal, 1:approximate cal, 2:accurate cal

# ---------- field solver.
# --- Setup r-z multigrid field solver 
#     Note: Grid and boundary conditions MUST be set before setup 
mesh_refine = false
w3d.solvergeom = w3d.RZgeom 
solver = MultiGrid2D() 
registersolver(solver)


# =================================================================================== #
# -----------------------------  set diagnostics  ----------------------------------- #
# =================================================================================== #

# set history diagnostic and moment accumulations 
top.nhist = 10                              # step interval for histories 
top.itmomnts[0:3] = [0,40,top.nhist]   # do loop ranges for moments 
                                            #   and status writes to tty

# set extrapolation parameter w to be zero
# to avoid bad extrapolation in moments calculation.
top.zmmntdtextmax = 0.

# invoke setup routine for graphics and output files (THIS IS MANDATORY)
setup()
 
# Plots of processed data 
# --- B_z(r,z) vs z at a few values of r 
plg(bz_m[0,      :],z_m/mm,color='black')
plg(bz_m[nr*0.25,:],z_m/mm,color='blue' )
plg(bz_m[nr*0.45,  :],z_m/mm,color='red'  )
limits(-450,450,0,1.1)
ptitles('Bz: r=0,R/2,09*R Black,Blue,Red','z [mm]','Bz(r,z)/Bz(0,0)')
fma()

# --- B_r(r,z) vs z at a few values of r
plg(br_m[1,:],z_m/mm,color='black')
plg(br_m[nr*0.25,:],z_m/mm,color='blue')
plg(br_m[nr*0.45,:],z_m/mm,color='red')
limits(-450,450,-0.5,0.5)
ptitles('Br: r=0,R/2,09*R Black,Blue,Red','z [mm]','Br(r,z)/Bz(0,0)')
fma()

# --- B_r(r,z) vs r at a few values of r
plg(bz_m[:,nz*0.5],r_m/mm,color='black')
plg(bz_m[:,nz*0.505],r_m/mm,color='blue')
plg(bz_m[:,nz*0.51],r_m/mm,color='red')
limits(0,90,0,1.1)
ptitles('Bz: z=0,Black,Blue,Red','r [mm]','Bz(r,z)/Bz(0,0)')
fma()

package("w3d"); generate()

# Uncomment to turn off space-charge deposition for single particle runs 
top.depos = "none"

# Install conducting rods on the mesh 
fieldsol()

# ---------- write out particle data
def diag_calls():
  if top.it%2==0:
    fdist = open("orbit1.dat", "a")
    fdist.write("%10d %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e\n"%(top.it, spi.getx()[1], spi.gety()[1], spi.getz()[1], spi.getvx()[1], spi.getvy()[1], spi.getvz()[1]))
    fdist.close()
    fdist = open("orbit2.dat", "a")
    fdist.write("%10d %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e\n"%(top.it, spi.getx()[2], spi.gety()[2], spi.getz()[2], spi.getvx()[2], spi.getvy()[2], spi.getvz()[2]))
    fdist.close()
    fdist = open("orbit3.dat", "a")
    fdist.write("%10d %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e\n"%(top.it, spi.getx()[3], spi.gety()[3], spi.getz()[3], spi.getvx()[3], spi.getvy()[3], spi.getvz()[3]))
    fdist.close()
    fdist = open("orbit0.dat", "a")
    fdist.write("%10d %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e\n"%(top.it, spi.getx()[0], spi.gety()[0], spi.getz()[0], spi.getvx()[0], spi.getvy()[0], spi.getvz()[0]))
    fdist.close()

# Install diagnostic calls after simulation step
installafterstep(diag_calls)

# Setup for timing 
wtimeon() 
step(nstep*nadvance)
fma()

# output data to a pdb file.
fout=PW("final.pdb")
fout.np = spi.getn()
fout.x = spi.getx()
fout.y = spi.gety()
fout.z = spi.getz()
fout.vx = spi.getvx()
fout.vy = spi.getvy()
fout.vz = spi.getvz()
fout.close()


# Print out timing statistics of run 
printtimers()
