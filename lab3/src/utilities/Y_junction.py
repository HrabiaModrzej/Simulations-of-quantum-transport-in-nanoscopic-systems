import kwant
import numpy as np
import matplotlib.pyplot as plt

#the calculations are done in atomic units e=h=me=1. Here the conversion factors are defined.
def eV2au(energy): #eV -> j.a
    return energy*0.03674932587122423
def au2eV(energy): #j.a -> eV
    return energy*27.2117
def nm2au(length): #nm -> j.a
    return length*18.89726133921252
def T2au(length):  #Tesla -> j.a
    return length*4.254382E-6

class Y_Junction():
    def __init__(self, *, dx = nm2au(2), L = int(100), W = int(15), m = 0.014, V0 = eV2au(0.05), sigma = nm2au(10),
                  B = T2au(0), R1 = nm2au(60), R2 = nm2au(120)):
        self.m  = m
        self.dx = dx
        self.L  = L
        self.W  = W
        self.V0 = V0
        self.B  = B
        self.R1 = R1
        self.R2 = R2
        self.sigma = sigma

        def is_y_junction_shape(pos):
            (x, y)  = pos
            dist_sq = x**2 + y**2

            if( x <= -0.5*(R2+R1) and x >= -R2-L*dx/2 and y >= -W*dx and y <= W*dx ):
                return True
            
            if( x >= 0 and x <= L*dx/2):
                if((y >= R1 and y <= R2) or (y >= -R2 and y <= -R1)):
                    return True
                
            if( x <= 0 and dist_sq >= R1**2 and dist_sq <= R2**2):
                return True
            
            return False
        
        self.is_the_shape = is_y_junction_shape

def make_system(yj):
    m  = yj.m
    dx = yj.dx
    L  = yj.L
    W  = yj.W
    V0 = yj.V0
    x0 = nm2au(L)
    y0 = 0.0
    B  = yj.B
    sigma = yj.sigma

    t = 1/(2*m*dx**2)

    def potential(x, y):
        return V0 * np.exp( ( - (x-x0)**2 - (y-y0)**2 ) / sigma**2 )

    def onsite(site):
        (x,y) = site.pos
        return 4*t + potential(x,y)

    def hopping(site_i, site_j):
        (x_i, y_i) = site_i.pos
        (x_j, y_j) = site_j.pos
        return -t * np.exp( -0.5j * B * (x_i - x_j) * (y_i + y_j) )
    
    lat = kwant.lattice.square(dx, norbs=1)
    sys = kwant.Builder()

    sys[lat.shape(yj.is_the_shape, (0, 0.5*(yj.R1 + yj.R2)))] = onsite
    sys[(kwant.builder.HoppingKind((-1,0), lat, lat))] = hopping
    sys[(kwant.builder.HoppingKind((0,-1), lat, lat))] = hopping

    lead_left = kwant.Builder(kwant.TranslationalSymmetry((-dx, 0)))    
    lead_left[(lat(-yj.L/2,j) for j in range(-W,W))] = onsite
    lead_left[(kwant.builder.HoppingKind((-1,0), lat, lat))] = hopping
    lead_left[(kwant.builder.HoppingKind((0,-1), lat, lat))] = hopping
    sys.attach_lead(lead_left)

    lead_top_right = kwant.Builder(kwant.TranslationalSymmetry((dx, 0)))    
    lead_top_right[(lat(yj.L/2,j + yj.R1/yj.dx) for j in range(2*W))] = onsite
    lead_top_right[(kwant.builder.HoppingKind((-1,0), lat, lat))] = hopping
    lead_top_right[(kwant.builder.HoppingKind((0,-1), lat, lat))] = hopping
    sys.attach_lead(lead_top_right)

    lead_bot_right = kwant.Builder(kwant.TranslationalSymmetry((dx, 0)))    
    lead_bot_right[(lat(yj.L/2,j - yj.R2/yj.dx) for j in range(2*W))] = onsite
    lead_bot_right[(kwant.builder.HoppingKind((-1,0), lat, lat))] = hopping
    lead_bot_right[(kwant.builder.HoppingKind((0,-1), lat, lat))] = hopping
    sys.attach_lead(lead_bot_right)

    sys = sys.finalized()
    return sys

def disperssion(yj, nr_lead, k_max, k_steps):
    dx  = yj.dx
    sys = make_system(yj)
    momenta = np.linspace(-k_max*dx, k_max*dx, k_steps)
    bands = kwant.physics.Bands(sys.leads[nr_lead])
    energies = [bands(k) for k in momenta]
    return (momenta/dx)*nm2au(1.0), energies

def transmission(yj, E, lead_in, lead_out):
    E   = eV2au(E)
    sys = make_system(yj)
    smatrix = kwant.smatrix(sys,E)
    t = smatrix.transmission(lead_in, lead_out)
    return t

def conductance(yj, Emax, E_steps, lead_in, lead_out):
    energies = np.linspace(0, Emax, E_steps)
    cond = [transmission(yj, E, lead_in, lead_out) for E in energies]
    return energies, cond

def current(yj, E, nr_lead, nr_mod, *, ax=None):
    E   = eV2au(E)
    sys = make_system(yj)
    current = kwant.operator.Current(sys).bind()
    psi = kwant.wave_function(sys, E)(nr_lead)
    curr = current(psi[nr_mod])
    kwant.plotter.current(sys, curr, ax=ax)
