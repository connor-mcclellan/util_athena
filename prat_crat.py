import astropy.constants as c
import astropy.units as u

# Input Parameters
basename = 'rk3_implicit_decoupled_isowind_1e11K_zeta-psi-anglefix_kff0'
ns_mass = 1.4 * u.M_sun.cgs
rhobase = 100 * u.g / u.cm**3
temp = 1e11 * u.K
redd = 1.1                      # Eddington ratio
rin = 1.2e6 * u.cm              # NS radius in cgs units
vinflow = 1e-5 # * c, flow velocity at inner boundary
mu = 1.0
dt = 1e-1

# Characteristic quantities
v0 = c.c.cgs  # np.sqrt(c.k_B.cgs * temp / mu / c.m_p.cgs)
temp0 = (v0**2 * mu * c.m_p.cgs / c.k_B.cgs).to(u.K)
crat = (c.c.cgs / v0).decompose()
arad = 4 * c.sigma_sb.cgs / c.c.cgs
kes = 0.2 * u.cm**2 * u.g**(-1)
rg = (2 * c.G.cgs * ns_mass / (c.c.cgs**2)).to(u.cm)
prat = (arad * temp0**4 * kes * rg / v0**2).decompose()
rho0 = (1 / kes / rg).to(u.g / u.cm**3)

athinput = (
'''<comment>
problem   = Flux-driven outflow from surface of NS
reference =
configure = --coord spherical_polar --prob=eddington_wind_dimensionless -hdf5 -implicit

<job>
'''
"problem_id = {} # problem ID: basename of output filenames".format(basename)+"\n"
'''
<output1>
id = prim
file_type  = tab      # Binary data dump
variable   = prim      # variables to be output
ghost_zones = true
'''
"dt       = {}".format(dt)+"\n"
'''
<output2>
id = uov
file_type = tab
variable = uov
ghost_zones = true
'''
"dt       = {}".format(dt)+"\n"
'''
<output3>
id = rad_mom
file_type  = tab      # Binary data dump
variable   = Fr0       # variables to be output
ghost_zones = true
'''
"dt       = {}".format(dt)+"\n"
'''
<output4>
id = rst
file_type  = rst
'''
"dt       = {}".format(dt)+"\n"
'''
<time>
cfl_number = 0.4                # The Courant, Friedrichs, & Lewy (CFL) Number
integrator = rk3
nlim       = -1                 # cycle limit
'''
"tlim       = {}                # time limit".format((50*rin / c.c.cgs / (rg/v0)).value)+""
'''
dt0        = 1.0e-9             # Initial time step

<mesh>
nx1        = 256                 # Number of zones in X1-direction
'''
"x1min      = {}                 # minimum value of X1 in gravitational radii".format((rin/rg).decompose())+"\n"
"x1max      = {}                 # maximum value of x1".format((10*rin/rg).decompose())+"\n"
'''
ix1_bc     = user               # inner-X1 boundary flag
ox1_bc     = outflow            # inner-X1 boundary flag
ix1_rad_bc = user
ox1_rad_bc = outflow
x1rat      = 1.05

nx2        = 1                  # Number of zones in X2-direction
x2min      = 0.0                # minimum value of X2
x2max      = 3.141592653589793  # maximum value of X2
ix2_bc     = periodic           # inner-X2 boundary flag
ox2_bc     = periodic           # inner-X2 boundary flag
ix2_rad_bc = periodic
ox2_rad_bc = periodic
nx3        = 1                  # Number of zones in X3-direction
x3min      = 0.0                # minimum value of X3
x3max      = 3.141592653589793  # maximum value of X3
ix3_bc     = periodic           # inner-X3 boundary flag
ox3_bc     = periodic           # inner-X3 boundary flag
ix3_rad_bc = periodic
ox3_rad_bc = periodic

<meshblock>
nx1         = 256
nx2         = 1
nx3         = 1

<hydro>
decouple_hydro = 1
gamma       = 1.6666666666667       # gamma = C_p/C_v
'''
"dfloor      = {}".format(1e-5* rhobase/rho0)+"\n"
"pfloor      = {}".format(1e-9 * rhobase/rho0)+"\n"
'''
\n<radiation>
n_frequency = 1
nmu         = 1 # dummy value, not used when angle_flag = 1
nzeta       = 16
npsi        = 0
angle_flag  = 1
'''
"Prat        = {}".format(prat)+"\n"
"Crat        = {}".format(crat)+"\n"
"\n"
"<problem>\n"
"rhobase  = {}".format(rhobase/rho0)+"\n"
"isotemp = {}".format(temp/temp0)+"\n"
"rg       = {}".format(rg.decompose().to(u.cm).value)+"\n"
"redd     = {}".format(redd)+"\n"
"vinflow  = {}".format(vinflow)+"\n"
)

with open("athinput.eddingtonwind", "w") as f:
    f.writelines(athinput)

print("rho0 = {}".format(rho0))
print("temp0 = {}".format(temp0))
print("x1min = {}".format(rin/rg))
print("\nFile 'athinput.eddingtonwind' saved.")
