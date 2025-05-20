import kwant
import matplotlib.pyplot as plt
import numpy as np
from utilities import spin_precessor as SP
from utilities import utils as u

plt.rcParams["font.family"] = "Latin Modern Math"

def Task1_0():
    sp  = SP.SpinPrecessor()
    sys = SP.make_system(sp)
    kwant.plot(sys, fig_size=(10,3), colorbar=False, show=False, num_lead_cells=2)
    plt.tight_layout()
    plt.savefig("./../plots/task1_system.pdf")
    plt.show()

def Task1_1():
    sp  = SP.SpinPrecessor()
    momenta, energies = SP.disperssion(sp, 0, .1, 500)
    plt.figure(figsize=(8,4))
    plt.plot(momenta, np.asarray(energies)/u.eV2au(1.0),'k-')
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.ylim((0,.2))  
    plt.xlim((-0.5,.5))  
    plt.xlabel("k [1/nm]",fontsize=15)
    plt.ylabel("E [eV]",fontsize=15)
    plt.tight_layout()
    plt.savefig("./../plots/task1_disperion.pdf")
    plt.show()

def Task1_2():

    # B = [1 0 0]
    sp  = SP.SpinPrecessor(Bx=u.T2au(1))
    momenta, energies = SP.disperssion(sp, 0, .1, 500)
    plt.figure(figsize=(6,4))
    plt.plot(momenta, np.asarray(energies)/u.eV2au(1.0),'k-')
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.ylim((0,.2))  
    plt.xlim((-0.5,.5))  
    plt.xlabel("k [1/nm]",fontsize=15)
    plt.ylabel("E [eV]",fontsize=15)
    plt.tight_layout()
    plt.savefig("./../plots/task1_disperion_x.pdf")
    plt.show()

    # B = [0 1 0]
    sp  = SP.SpinPrecessor(By=u.T2au(1))
    momenta, energies = SP.disperssion(sp, 0, .1, 500)
    plt.figure(figsize=(6,4))
    plt.plot(momenta, np.asarray(energies)/u.eV2au(1.0),'k-')
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.ylim((0,.2))  
    plt.xlim((-0.5,.5))  
    plt.xlabel("k [1/nm]",fontsize=15)
    plt.ylabel("E [eV]",fontsize=15)
    plt.tight_layout()
    plt.savefig("./../plots/task1_disperion_y.pdf")
    plt.show()

    # B = [0 0 1]
    sp  = SP.SpinPrecessor(Bz=u.T2au(1))
    momenta, energies = SP.disperssion(sp, 0, .1, 500)
    plt.figure(figsize=(8,4))
    plt.plot(momenta, np.asarray(energies)/u.eV2au(1.0),'k-')
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.ylim((0,.2))  
    plt.xlim((-0.5,.5))  
    plt.xlabel("k [1/nm]",fontsize=15)
    plt.ylabel("E [eV]",fontsize=15)
    plt.tight_layout()
    plt.savefig("./../plots/task1_disperion_z.pdf")
    plt.show()

def Task1_3():
    sp  = SP.SpinPrecessor(Bz=u.T2au(1))
    energies, conductances = SP.conductance(sp, 0.05, 100)
    plt.figure(figsize=(8,4))
    plt.plot(energies, conductances,'k-')
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.xlabel("E (eV)",fontsize=15)
    plt.ylabel("G (e^2/h)",fontsize=15)
    plt.tight_layout()
    plt.savefig("./../plots/task1_conductance.pdf")
    plt.show()

def Task1_4():
    E = 5e-3
    By_vec = np.linspace(0,1,100)
    transmissions_00 = np.zeros(len(By_vec))
    transmissions_01 = np.zeros(len(By_vec))
    transmissions_10 = np.zeros(len(By_vec))
    transmissions_11 = np.zeros(len(By_vec))
    for i in range(len(By_vec)):
        sp  = SP.SpinPrecessor(Bz=u.T2au(0.1), add_By=u.T2au(By_vec[i]))
        transmissions_00[i] = SP.transmission_spin(sp, E, 0, 0)
        transmissions_01[i] = SP.transmission_spin(sp, E, 0, 1)
        transmissions_10[i] = SP.transmission_spin(sp, E, 1, 0)
        transmissions_11[i] = SP.transmission_spin(sp, E, 1, 1)
    
    plt.figure(figsize=(8,4))
    plt.plot(By_vec, transmissions_00, label='T_{up→up}')
    plt.plot(By_vec, transmissions_01, label='T_{up→dn}')
    plt.plot(By_vec, transmissions_10, label='T_{dn→up}', linestyle="--")
    plt.plot(By_vec, transmissions_11, label='T_{dn→dn}', linestyle="--")
    plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand", ncol=4)
    plt.tick_params(axis='both', which='major', labelsize=15)
    #plt.ylim((-0.5, 8.5))  
    #plt.xlim((-0.005, 0.055))  
    plt.xlabel("By (T)",fontsize=15)
    plt.ylabel("T",fontsize=15)
    plt.tight_layout()
    plt.savefig("./../plots/task1_transmissions_spin_dep.pdf")
    plt.show()

def Task1_5():
    E = u.eV2au(5e-3)
    sp  = SP.SpinPrecessor(Bz=u.T2au(0.1), add_By=u.T2au(0.6))
    fig, subfigs = plt.subplots(2, 1, figsize=(10, 4), dpi = 100)
    SP.density_up_map(sp, E, ax = subfigs[0])
    subfigs[0].set_title("Spin up", fontsize=15)
    SP.density_dn_map(sp, E, ax = subfigs[1])
    subfigs[1].set_title("Spin down", fontsize=15)
    plt.tight_layout()
    plt.savefig("../plots/task1_charge_density.pdf")
    plt.show()

def Task1_6():
    E = u.eV2au(5e-3)
    sp  = SP.SpinPrecessor(Bz=u.T2au(0.1), add_By=u.T2au(0.6))
    fig, subfigs = plt.subplots(3, 1, figsize=(10, 4), dpi = 100)
    SP.spin_x_density(sp, E, ax = subfigs[0])
    subfigs[0].set_title("s_x", fontsize=15)
    SP.spin_y_density(sp, E, ax = subfigs[1])
    subfigs[1].set_title("s_y", fontsize=15)
    SP.spin_z_density(sp, E, ax = subfigs[2])
    subfigs[2].set_title("s_z", fontsize=15)
    plt.tight_layout()
    plt.savefig("../plots/task1_spin_density.pdf")
    plt.show()

def Task2_1():
    alpha = u.nm2au(u.eV2au(18e-3))
    sp  = SP.SpinPrecessor(L=u.nm2au(800), W=u.nm2au(100), alpha=alpha)
    momenta, energies = SP.disperssion(sp, 0, .1, 500)
    plt.figure(figsize=(8,4))
    plt.plot(momenta, np.asarray(energies)/u.eV2au(1.0),'-k')
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.ylim((0,.05))  
    plt.xlim((-0.2,.2))  
    plt.xlabel("k [1/nm]",fontsize=15)
    plt.ylabel("E [eV]",fontsize=15)
    plt.tight_layout()
    plt.savefig("./../plots/task2_disperion.pdf")
    plt.show()

def Task2_2():
    alpha = u.nm2au(u.eV2au(18e-3))
    sp  = SP.SpinPrecessor(L=u.nm2au(800), W=u.nm2au(100), alpha=alpha)
    energies, conductances = SP.conductance(sp, 0.05, 100)
    plt.figure(figsize=(8,4))
    plt.plot(energies, conductances,'k-')
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.xlabel("E (eV)",fontsize=15)
    plt.ylabel("G (e^2/h)",fontsize=15)
    plt.tight_layout()
    plt.savefig("./../plots/task2_conductance.pdf")
    plt.show()

def Task2_3():
    E = 5e-3
    alpha_vec = np.linspace(0, 50, 100)
    transmissions_00 = np.zeros(len(alpha_vec))
    transmissions_01 = np.zeros(len(alpha_vec))
    transmissions_10 = np.zeros(len(alpha_vec))
    transmissions_11 = np.zeros(len(alpha_vec))
    for i in range(len(alpha_vec)):
        sp  = SP.SpinPrecessor(L=u.nm2au(800), W=u.nm2au(100), add_so=u.eV2au(u.nm2au(alpha_vec[i]*1e-3)))
        transmissions_00[i] = SP.transmission_spin(sp, E, 0, 0)
        transmissions_01[i] = SP.transmission_spin(sp, E, 0, 1)
        transmissions_10[i] = SP.transmission_spin(sp, E, 1, 0)
        transmissions_11[i] = SP.transmission_spin(sp, E, 1, 1)
    
    plt.figure(figsize=(8,4))
    plt.plot(alpha_vec, transmissions_00, label='T_{up→up}')
    plt.plot(alpha_vec, transmissions_01, label='T_{up→dn}')
    plt.plot(alpha_vec, transmissions_10, label='T_{dn→up}', linestyle="--")
    plt.plot(alpha_vec, transmissions_11, label='T_{dn→dn}', linestyle="--")
    plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand", ncol=4)
    plt.tick_params(axis='both', which='major', labelsize=15)
    #plt.ylim((-0.5, 8.5))  
    #plt.xlim((-0.005, 0.055))  
    plt.xlabel("alpha (meVnm)",fontsize=15)
    plt.ylabel("T",fontsize=15)
    plt.tight_layout()
    plt.savefig("./../plots/task2_transmissions_spin_dep.pdf")
    plt.show()

def Task2_4():
    P_vec = [0.2, 0.4, 1.0]
    alpha_vec = np.linspace(0, 50, 100)
    cond_up_vec = np.zeros(len(alpha_vec))
    cond_dn_vec = np.zeros(len(alpha_vec))
    cond_sum_vec= np.zeros(len(alpha_vec))
    fig, subfigs = plt.subplots(1, 3, figsize=(15, 4), dpi = 100)
    for j in range(len(P_vec)):
        print("Rozpoczynam obliczenia dla P="+str(P_vec[j]))
        for i in range(len(alpha_vec)):
            sp  = SP.SpinPrecessor(L=u.nm2au(800), W=u.nm2au(100), add_so=u.eV2au(u.nm2au(alpha_vec[i]*1e-3)))
            cond_dn_vec[i] = SP.spin_conductance_dn(sp,P_vec[j])
            cond_up_vec[i] = SP.spin_conductance_up(sp,P_vec[j])
            cond_sum_vec[i] = SP.spin_conductance_sum(P_vec[j], cond_up_vec[i], cond_dn_vec[i])
        subfigs[j].plot(alpha_vec, cond_up_vec, label='G_{up}')
        subfigs[j].plot(alpha_vec, cond_dn_vec, label='G_{dn}')
        subfigs[j].plot(alpha_vec, cond_sum_vec, label='G')
        subfigs[j].legend(bbox_to_anchor=(0,-0.4,1,0.2), loc="lower left", mode="expand", ncol=3)
        subfigs[j].tick_params(axis='both', which='major', labelsize=12)
        subfigs[j].set_xlabel("alpha (meVnm)",fontsize=12)
        subfigs[j].set_ylabel("G",fontsize=12)
        subfigs[j].set_title("P="+str(P_vec[j]), fontsize=15)
        print("Kończę obliczenia dla P="+str(P_vec[j]))
    plt.tight_layout()
    plt.savefig("./../plots/task2_conds_spin_dep.pdf")
    plt.show()

def Task2_5():
    E = u.eV2au(5e-3)
    alpha = u.nm2au(u.eV2au(18e-3))

    sp  = SP.SpinPrecessor(L=u.nm2au(800), W=u.nm2au(100), add_so=alpha)
    fig, subfigs = plt.subplots(2, 1, figsize=(10, 4), dpi = 100)
    SP.density_up_map(sp, E, ax = subfigs[0])
    subfigs[0].set_title("Spin up", fontsize=15)
    SP.density_dn_map(sp, E, ax = subfigs[1])
    subfigs[1].set_title("Spin down", fontsize=15)
    plt.tight_layout()
    plt.savefig("../plots/task2_charge_density.pdf")
    plt.show()

    fig, subfigs = plt.subplots(3, 1, figsize=(10, 4), dpi = 100)
    SP.spin_x_density(sp, E, ax = subfigs[0])
    subfigs[0].set_title("s_x", fontsize=15)
    SP.spin_y_density(sp, E, ax = subfigs[1])
    subfigs[1].set_title("s_y", fontsize=15)
    SP.spin_z_density(sp, E, ax = subfigs[2])
    subfigs[2].set_title("s_z", fontsize=15)
    plt.tight_layout()
    plt.savefig("../plots/task2_spin_density.pdf")
    plt.show()


Task2_2()

