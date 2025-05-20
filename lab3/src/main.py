import kwant
import numpy as np
import matplotlib.pyplot as plt
from utilities import nanowire as NW
from utilities import Y_junction as YJ

plt.rcParams["font.family"] = "Latin Modern Math"


def Task1():
    nw = NW.NanoWire()
    sys = NW.make_system(nw)
    kwant.plot(sys, site_color=lambda site: sys.hamiltonian(site,site), fig_size=(10,5), colorbar=False, show=False, num_lead_cells=2)
    plt.savefig("../plots/task1_system.pdf")
    #plt.show()
    print("System defined")

    print("Starting the conductance calculations")
    energies, conductances = NW.conductance(nw, 0.2, 50)
    plt.figure()
    plt.plot(energies, conductances, color = 'black')
    plt.xlabel("E [eV]")
    plt.ylabel("G [2e^2/h]")
    plt.savefig("../plots/task1_condutance.pdf")
    #plt.show()

    print("Plotting wavefunctions and currents")
    fig, subfigs = plt.subplots(2, 3, figsize=(20, 5), dpi = 100)
    NW.wave_function(nw, 0.03, 0, ax = subfigs[0, 0])
    subfigs[0, 0].set_title("Energy = 0.03 eV", fontsize=15)
    NW.wave_function(nw, 0.05, 0, ax = subfigs[0, 1])
    subfigs[0, 1].set_title("Energy = 0.05 eV", fontsize=15)
    NW.wave_function(nw, 0.10, 0, ax = subfigs[0, 2])
    subfigs[0, 2].set_title("Energy = 0.1 eV",  fontsize=15)
    NW.current(nw, 0.03, 0, 0, ax = subfigs[1, 0])
    NW.current(nw, 0.05, 0, 0, ax = subfigs[1, 1])
    NW.current(nw, 0.10, 0, 0, ax = subfigs[1, 2])
    plt.tight_layout()
    plt.savefig("../plots/task1_wavefunctions_currents.pdf")
    plt.show()

def Task2():
    nw  = NW.NanoWire(V0 = 0.0, B = NW.T2au(2.0), W = int(20))
    sys = NW.make_system(nw)
    kwant.plot(sys, site_color=lambda site: sys.hamiltonian(site,site), fig_size=(10,5), colorbar=False, show=False, num_lead_cells=2)
    plt.savefig("../plots/task2_system.pdf")
    #plt.show()

    moments, enes = NW.disperssion(nw, 0, .1, 200)
    plt.figure(dpi=100)
    plt.plot(moments, np.asarray(enes)/NW.eV2au(1.0),'k-')
    plt.ylim((0,.2))  
    plt.xlim((-0.5,.5))  
    plt.xlabel("k [1/nm]",fontsize=12)
    plt.ylabel("E [eV]",fontsize=12)
    plt.tight_layout()
    plt.savefig("../plots/task2_disp_rel_B2T_W40.pdf")
    #plt.show()
    
    nw  = NW.NanoWire(V0 = 0.0, B = NW.T2au(2.0), W = int(50))
    moments, enes = NW.disperssion(nw, 0, .1, 200)
    plt.figure(dpi=100)
    plt.plot(moments, np.asarray(enes)/NW.eV2au(1.0),'k-')
    plt.ylim((0,.2))  
    plt.xlim((-0.5,.5))  
    plt.xlabel("k [1/nm]",fontsize=12)
    plt.ylabel("E [eV]",fontsize=12)
    plt.tight_layout()
    plt.savefig("../plots/task2_disp_rel_B2T_W100.pdf")
    #plt.show()

    energies, conductances = NW.conductance(nw, 0.2, 100)
    plt.figure()
    plt.plot(energies, conductances, color = 'black')
    plt.xlabel("E [eV]")
    plt.ylabel("G [2e^2/h]")
    plt.savefig("../plots/task2_condutance.pdf")
    #plt.show()

    fig, ax = plt.subplots(1, 1, dpi=100)
    NW.wave_function(nw, 0.015, 0, ax=ax)
    ax.set_xlabel("x [a.u.]")
    ax.set_ylabel("y [a.u.]")
    plt.savefig("../plots/task2_wavefunction_left.pdf")
    plt.show()
    
    fig, ax = plt.subplots(1, 1, dpi=100)
    NW.wave_function(nw, 0.015, 1, ax=ax)
    ax.set_xlabel("x [a.u.]")
    ax.set_ylabel("y [a.u.]")
    plt.savefig("../plots/task2_wavefunction_right.pdf")
    plt.show()

def Task3():
    yj = YJ.Y_Junction(V0 = 0., B = YJ.T2au(0), W = int(15), L = int(50))
    sys = YJ.make_system(yj)
    kwant.plot(sys, site_color=lambda site: sys.hamiltonian(site,site), fig_size=(10,5), colorbar=False, show=False, num_lead_cells=2);
    plt.savefig("../plots/task3_system.pdf")
    #plt.show()

    k_points, energies = YJ.disperssion(yj, 0, 0.1, 200)
    plt.figure(figsize = (10,7))
    plt.plot(k_points, np.asarray(energies)/YJ.eV2au(1.0),'k-')
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.ylim((0,.2))  
    plt.xlim((-0.5,.5))  
    plt.xlabel("k [1/nm]",fontsize=15)
    plt.ylabel("E [eV]",fontsize=15)
    plt.savefig("../plots/task3_dispersion_relation.pdf")
    #plt.show()

    B_vec = np.linspace(YJ.T2au(-10), YJ.T2au(10), 50)
    conds_up_lead = np.zeros(len(B_vec))
    conds_dn_lead = np.zeros(len(B_vec))
    E = 0.1

    for i in range(len(B_vec)):
        yj = YJ.Y_Junction(V0 = 0, B = B_vec[i], W = 15, L = 50)
        conds_up_lead[i] = YJ.transmission(yj, E, 1, 0)
        conds_dn_lead[i] = YJ.transmission(yj, E, 2, 0)

    plt.figure()
    plt.plot(B_vec, conds_up_lead, color = '#00693C', label = "Conductance upper lead")
    plt.plot(B_vec, conds_dn_lead, color = '#A71930', label = "Conductance bottom lead")
    plt.xlabel("Bz [T]")
    plt.ylabel("G [2e^2/h]")
    plt.legend()
    plt.tight_layout()
    plt.savefig("../plots/task3_conductances.pdf")
    #plt.show()

    figs, subfigs = plt.subplots(2, 2, dpi=100)
    
    B = 0
    yj = YJ.Y_Junction(V0 = 0, B = YJ.T2au(B), W = int(15), L = int(50))
    YJ.current(yj, 0.05, 0, 0, ax = subfigs[0,0])
    subfigs[0,0].set_xlabel("x [a.u.]")
    subfigs[0,0].set_ylabel("y [a.u.]")
    subfigs[0,0].set_title(f"B = {B} T")

    # smooth conductance
    B = 1
    yj = YJ.Y_Junction(V0 = 0, B = YJ.T2au(B), W = int(15), L = int(50))
    YJ.current(yj, 0.05, 0, 0, ax = subfigs[0,1])
    subfigs[0,1].set_xlabel("x [a.u.]")
    subfigs[0,1].set_ylabel("y [a.u.]")
    subfigs[0,1].set_title(f"B = {B} T")

    # plateau
    B = 3
    yj = YJ.Y_Junction(V0 = 0, B = YJ.T2au(B), W = int(15), L = int(50))
    YJ.current(yj, 0.05, 0, 0, ax = subfigs[1,0])
    subfigs[1,0].set_xlabel("x [a.u.]")
    subfigs[1,0].set_ylabel("y [a.u.]")
    subfigs[1,0].set_title(f"B = {B} T")
    
    B = 5
    yj = YJ.Y_Junction(V0 = 0, B = YJ.T2au(B), W = int(15), L = int(50))
    YJ.current(yj, 0.05, 0, 0, ax = subfigs[1,1])
    subfigs[1,1].set_xlabel("x [a.u.]")
    subfigs[1,1].set_ylabel("y [a.u.]")
    subfigs[1,1].set_title(f"B = {B} T")

    plt.tight_layout()
    plt.savefig("../plots/task3_current.pdf")
    #plt.show()


#Task1()
#Task2()
Task3()
