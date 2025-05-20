import kwant
import numpy as np
from . import utils as u

mu_B = 0.5
g = -50

sigma_x = np.array([[0,1],[1,0]])
sigma_y = np.array([[0,-1j],[1j,0]])
sigma_z = np.array([[1,0],[0,-1]])
identity = np.array([[1,0],[0,1]])

class SpinPrecessor():
    def __init__(self, *, dx = u.nm2au(4), L = u.nm2au(2000), W = u.nm2au(100), m = 0.014, \
                 Bx = u.T2au(0), By = u.T2au(0), Bz = u.T2au(0), alpha = u.nm2au(0), add_By = u.T2au(0), add_so = 0):
        self.dx = dx
        self.L  = int(L/dx)
        self.W  = int(W/dx)
        self.m  = m
        self.Bx = Bx
        self.By = By
        self.Bz = Bz
        self.alpha = alpha
        self.add_By = add_By
        self.add_so = add_so

def make_system(ns):
    m     = ns.m
    dx    = ns.dx
    alpha = ns.alpha
    Bx    = ns.Bx
    By    = ns.By
    Bz    = ns.Bz
    L     = ns.L
    W     = ns.W
    add_By= ns.add_By
    add_so= ns.add_so

    t = 1/(2*m*dx**2)
    t_so = alpha/(2*dx)
    t_add = add_so/(2*dx)

    def onsite(site):
        (x,y) = site.pos
        if(x >= 0.2*L*dx and x<= 0.8*L*dx):
            return 4*t*identity + 0.5*mu_B*g*(Bx*sigma_x + add_By*sigma_y + By*sigma_y + Bz*sigma_z)
        return 4*t*identity + 0.5*mu_B*g*(Bx*sigma_x + By*sigma_y + Bz*sigma_z)
    
    def hopping_x(sitei, sitej):
        (xi, yi) = sitei.pos
        if (xi >= 0.2*L*dx and xi<= 0.8*L*dx):
            return -1*t*identity + 1j*t_so*sigma_y + 1j*t_add*sigma_y
        return -1*t*identity + 1j*t_so*sigma_y 
    
    def hopping_y(sitei, sitej):
        (xi, yi) = sitei.pos
        if (xi >= 0.2*L*dx and xi<= 0.8*L*dx):
            return -1*t*identity + 1j*t_so*sigma_x + 1j*t_add*sigma_x
        return -1*t*identity - 1j*t_so*sigma_x 

    lat = kwant.lattice.square(dx, norbs=2)
    sys = kwant.Builder()
    sys[(lat(i,j) for i in range(L) for j in range(W))]=onsite
    sys[(kwant.builder.HoppingKind((-1,0), lat, lat))] = hopping_x
    sys[(kwant.builder.HoppingKind((0,-1), lat, lat))] = hopping_y

    # I've changed it and the results are the same
    sigma_law = np.array([[-1,0],[0,1]])

    # Conservation law in both leads → crucial for spin dependent transmission calc
    left_lead = kwant.Builder(kwant.TranslationalSymmetry((-dx, 0)),conservation_law=sigma_law)
    left_lead[(lat(0,j) for j in range(W))]=onsite
    left_lead[(kwant.builder.HoppingKind((-1,0), lat, lat))] = hopping_x
    left_lead[(kwant.builder.HoppingKind((0,-1), lat, lat))] = hopping_y
    sys.attach_lead(left_lead)

    right_lead = kwant.Builder(kwant.TranslationalSymmetry((dx, 0)),conservation_law=sigma_law)    
    right_lead[(lat(0,j) for j in range(W))]=onsite
    right_lead[(kwant.builder.HoppingKind((-1,0), lat, lat))] = hopping_x
    right_lead[(kwant.builder.HoppingKind((0,-1), lat, lat))] = hopping_y
    sys.attach_lead(right_lead)
    
    #finalize the system
    sys = sys.finalized()
    return sys

#calculates the dispersion relation in the contact nr_lead in the range [-k_max,k_max] with nk points
def disperssion(nw, nr_lead, k_max, nk):
    dx=nw.dx
    sys=make_system(nw)
    momenta = np.linspace(-k_max*dx,k_max*dx,nk)
    bands=kwant.physics.Bands(sys.leads[nr_lead])
    energies=[bands(k) for k in momenta]
    return (momenta/dx)*u.nm2au(1.0), energies

#calculates the transmission coefficient
def transmission(nw, E):
    E=u.eV2au(E)
    sys=make_system(nw)
    smatrix=kwant.smatrix(sys,E)
    t=smatrix.transmission(1,0)
    return t

def transmission_spin(nw, E, spinin:int, spinout:int):
    E=u.eV2au(E)
    sys=make_system(nw)
    smatrix=kwant.smatrix(sys,E)
    t=smatrix.transmission((1, spinout), (0, spinin))
    return t

#calculates the conductance - the Landauer formula is used
def conductance(nw, Emax, ne):
    energies=np.linspace(0,Emax,ne)
    cond=[transmission(nw, E) for E in energies]
    return energies, cond

def density_up_map(NanoSystem, E, *, nr_lead=0, nr_mode =0 , ax=None):
    sys=make_system(NanoSystem)
    density_up = np.array([[1,0],[0,0]])
    density_up_op = kwant.operator.Density(sys,density_up)
    psi=kwant.wave_function(sys,E)(nr_lead)
    density_map = density_up_op(psi[nr_mode])
    kwant.plotter.map(sys, density_map, ax=ax)

def density_dn_map(NanoSystem, E, *, nr_lead=0, nr_mode =0 , ax=None):
    sys=make_system(NanoSystem)
    density_dn=np.array([[0,0],[0,1]])
    density_dn_op = kwant.operator.Density(sys,density_dn)
    psi=kwant.wave_function(sys,E)(nr_lead)
    density_map = density_dn_op(psi[nr_mode])
    kwant.plotter.map(sys, density_map, ax=ax)

def density_both_map(NanoSystem, E, *, nr_lead=0, nr_mode =0 , ax=None):
    sys=make_system(NanoSystem)
    density_both=np.array([[1,0],[0,1]])
    density_both_op = kwant.operator.Density(sys,density_both)
    psi=kwant.wave_function(sys,E)(nr_lead)
    density_map = density_both_op(psi[nr_mode])
    kwant.plotter.map(sys, density_map, ax=ax)

def spin_x_density(NanoSystem, E, *,nr_lead=0, nr_mode =0 , ax=None):
    sys=make_system(NanoSystem)
    spin_den_op = kwant.operator.Density(sys, sigma_x)
    psi=kwant.wave_function(sys,E)(nr_lead)
    density_map = spin_den_op(psi[nr_mode])
    kwant.plotter.map(sys, density_map, cmap='RdBu',ax=ax)

def spin_y_density(NanoSystem, E, *,nr_lead=0, nr_mode =0 , ax=None):
    sys=make_system(NanoSystem)
    spin_den_op = kwant.operator.Density(sys, sigma_y)
    psi=kwant.wave_function(sys,E)(nr_lead)
    density_map = spin_den_op(psi[nr_mode])
    kwant.plotter.map(sys, density_map, cmap='RdBu',ax=ax)

def spin_z_density(NanoSystem, E, *,nr_lead=0, nr_mode =0 , ax=None):
    sys=make_system(NanoSystem)
    spin_den_op = kwant.operator.Density(sys, sigma_z)
    psi=kwant.wave_function(sys,E)(nr_lead)
    density_map = spin_den_op(psi[nr_mode])
    kwant.plotter.map(sys, density_map, cmap='RdBu', ax=ax)

def spin_conductance_up(NanoSystem, P):
    E = 5e-3
    return 0.5*(1+P)*transmission_spin(NanoSystem,E,0,0) +\
    0.5*(1-P)*transmission_spin(NanoSystem,E,1,0)

def spin_conductance_dn(NanoSystem, P):
    E = 5e-3
    return 0.5*(1+P)*transmission_spin(NanoSystem,E,0,1) +\
    0.5*(1-P)*transmission_spin(NanoSystem,E,1,1)

def spin_conductance_sum(P, Gup, Gdn):
    return 0.5*(1+P)*Gup + 0.5*(1-P)*Gdn