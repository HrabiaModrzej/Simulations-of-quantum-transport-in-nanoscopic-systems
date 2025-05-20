import kwant 
import numpy as np
import matplotlib.pyplot as plt

#conversion functions into atomic units
def eV2au(energy): #eV -> j.a
    return energy*0.03674932587122423
def au2eV(energy): #j.a -> eV
    return energy*27.2117
def nm2au(length): #nm -> j.a
    return length*18.89726133921252
def T2au(length):  #T -> j.a
    return length*4.254382E-6


def make_system(nw):
    # Values used for defining the system
    m  = nw.m
    dx = nw.dx
    L  = nw.L
    W  = nw.W
    V0 = nw.V0
    x0 = nm2au(L)
    y0 = 0
    B  = nw.B 
    sigma = nw.sigma

    # Value showing up in Hamiltonian t = ℏ/(2mΔx²)
    t = 1/(2*m*dx**2)

    # The potential function in a system
    def potential(x, y):
        return V0 * np.exp( ( - (x-x0)**2 - (y-y0)**2 ) / sigma**2 )
    
    # Diagonal elements in defined finite differences 
    def onsite(site):
        (x,y) = site.pos
        return 4*t + potential(x,y)
    
    # Non-diagonal elements → hopping refers to the changes of the energy when moving the particle to the different site
    def hopping(site_i, site_j):
        (x_i, y_i) = site_i.pos
        (x_j, y_j) = site_j.pos
        return -t * np.exp( -0.5j * B * (x_i - x_j) * (y_i + y_j) )
    
    # Defining a system
    sys = kwant.Builder()
    # Describing the grid
    lat = kwant.lattice.square(dx, norbs=1)

    # Defining the hamiltonian for system
    sys[( lat(i,j) for i in range(L) for j in range(-W,W) )] = onsite
    sys[(kwant.builder.HoppingKind( (-1,0), lat, lat ))] = hopping
    sys[(kwant.builder.HoppingKind( (0,-1), lat, lat ))] = hopping

    # Defining and attaching the left contact(lead) to the system
    left_lead = kwant.Builder(kwant.TranslationalSymmetry((-dx,0)))
    left_lead[( lat(0,j) for j in range(-W,W) )] = onsite
    left_lead[(kwant.builder.HoppingKind( (-1,0), lat, lat ))] = hopping
    left_lead[(kwant.builder.HoppingKind( (0,-1), lat, lat ))] = hopping
    sys.attach_lead(left_lead)

    # Same for the right one
    right_lead = kwant.Builder(kwant.TranslationalSymmetry((dx,0))) # change -dx → dx translation out of the system
    right_lead[( lat(0,j) for j in range(-W,W) )] = onsite
    right_lead[(kwant.builder.HoppingKind( (-1,0), lat, lat ))] = hopping
    right_lead[(kwant.builder.HoppingKind( (0,-1), lat, lat ))] = hopping
    sys.attach_lead(right_lead)

    # Finalization 
    sys = sys.finalized()
    return sys


#calculates the dispersion relation in the contact nr_lead in the range [-k_max,k_max] with nk points
def disperssion(nw, nr_lead, k_max, k_steps):
    dx  = nw.dx
    sys = make_system(nw)
    momenta = np.linspace(-k_max*dx, k_max*dx, k_steps)
    bands = kwant.physics.Bands(sys.leads[nr_lead])
    energies = [bands(k) for k in momenta]
    return (momenta/dx)*nm2au(1.0), energies

#calculates the reflection and transmission coefficient
def transmission_reflection(nw, E):
    E = eV2au(E)
    sys = make_system(nw)
    smatrix = kwant.smatrix(sys, E)
    r = smatrix.transmission(0,0)
    t = smatrix.transmission(1,0)
    return r, t

#calculates the transmission coefficient
def transmission(nw, E):
    E   = eV2au(E)
    sys = make_system(nw)
    smatrix = kwant.smatrix(sys,E)
    t = smatrix.transmission(1,0)
    return t


#calculates the conductance - the Landauer formula is used
def conductance(nw, Emax, E_steps):
    energies = np.linspace(0, Emax, E_steps)
    cond = [transmission(nw, E) for E in energies]
    return energies, cond


#plots the wave function of an electron with energy E incident in the contact nr_lead
def wave_function(nw, E, nr_lead, *, ax= None):
    E    = eV2au(E)
    sys  = make_system(nw)
    wave = kwant.wave_function(sys, E)
    density = (abs(wave(nr_lead))**2).sum(axis=0)
    kwant.plotter.map(sys, density, ax=ax)

#plots the dos of an electron with energy E
def dos(nw, E):
    E   = eV2au(E)
    sys = make_system(nw)
    dos = kwant.ldos(sys, E)
    f   = kwant.plotter.map(sys,dos)
    return f

#plots the current of an electron with energy E incident in the contact nr_lead in the state nr_mod
def current(nw, E, nr_lead, nr_mod, *, ax=None):
    E   = eV2au(E)
    sys = make_system(nw)
    current = kwant.operator.Current(sys).bind()
    psi = kwant.wave_function(sys, E)(nr_lead)
    curr = current(psi[nr_mod])
    kwant.plotter.current(sys, curr, ax=ax)


class NanoWire():
    def __init__(self, *, dx = nm2au(2), L = int(100), W = int(15), m = 0.014
                        , V0 = eV2au(0.05), sigma = nm2au(10), B = T2au(0)):
        self.dx = dx
        self.L  = int(L)
        self.W  = int(W)
        self.m  = m
        self.V0 = V0
        self.B  = B
        self.sigma = sigma