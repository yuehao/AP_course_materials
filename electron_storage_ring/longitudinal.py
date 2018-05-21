
import numpy as np
import math


def longitudinal_evolve_delta(turns, phi_list_ini, delta_list_ini, sin_phi_s=0, E0_ini=100e9, mass=938e6, e_volt=5e6, alphac=0.002, harm=360.0, E_change=False, update_eta=True):
    gamma0 = E0_ini / mass
    P0_old = math.sqrt(E0_ini * E0_ini - mass * mass)
    beta0 = math.sqrt(1 - 1.0 / gamma0 / gamma0)
    phi_list=[]
    phi_list.append(phi_list_ini)
    delta_list=[]
    delta_list.append(delta_list_ini)
    E0=E0_ini
    #nus = math.sqrt(harm * abs(ita) * e_volt / 2 / np.pi / E0_ini / beta / beta)
    for _ in range(turns):
        #yield (phi_list, delta_list)
        pl=phi_list[-1]*1.0
        dl=delta_list[-1]*1.0
        dl +=  e_volt * (np.sin(pl) - sin_phi_s ) / beta0 /beta0 / E0
        if update_eta:
            delta_gamma = dl * beta0 * beta0 * E0 / mass
        else:
            delta_gamma =0
        eta = alphac - 1.0/(gamma0+delta_gamma)/(gamma0+delta_gamma)
        pl += 2.0*np.pi*harm*eta*dl

        if E_change:

            E0 += e_volt*sin_phi_s
            P0 = math.sqrt(E0 * E0 - mass * mass)
            gamma0 = E0 / mass
            beta0 = math.sqrt(1 - 1.0 / gamma0 / gamma0)
            dl*=(P0_old/P0)
            P0_old=P0
        phi_list.append(pl)
        delta_list.append(dl)

    return np.vstack(phi_list), np.vstack(delta_list)



def longitudinal_evolve(turns, phi_list_ini, dE_list_ini, sin_phi_s=0, E0_ini=100e9, mass=938e6, e_volt=5e6, alphac=0.002, harm=360.0, update_eta=True, energy_change=False,
                        damping_turn=0, excitation=False, sr_power=0):
    E0 = E0_ini
    p0 = np.sqrt(E0 * E0 - mass * mass)
    p0_ini=np.sqrt(E0_ini*E0_ini-mass*mass)
    gamma0 = E0 / mass
    beta0 = math.sqrt(1 - 1.0 / gamma0 / gamma0)
    phi_list=[]
    phi_list.append(phi_list_ini)

    dE_list=[]
    dE_list.append(dE_list_ini)

    delta_list=[]
    e_temp=np.array(dE_list_ini)+E0

    dl_ini= np.sqrt(e_temp*e_temp-mass*mass)/p0_ini-1.0
    delta_list.append(dl_ini)

    #nus = math.sqrt(harm * abs(ita) * e_volt / 2 / np.pi / E0_ini / beta / beta)

    for _ in range(turns):
        #yield (phi_list, delta_list)
        pl=phi_list[-1]*1.0
        dEl=dE_list[-1]*1.0
        dEl +=  e_volt * (np.sin(pl) - sin_phi_s )
        if excitation:
            alpha=0.0072973525664
            nph=5*np.pi*alpha*gamma0
            srp=sr_power/np.power(E0, 4.0)*np.power(E0+dEl, 4.0)
            ave_u=srp/nph
            ave_u0=sr_power/nph


        if damping_turn>0:
            dEl *= (1-1.0/damping_turn)
        if energy_change:
            E0 += e_volt * sin_phi_s
            p0 = np.sqrt(E0 * E0 - mass * mass)
        dl = np.sqrt((E0 + dEl) * (E0 + dEl) - mass * mass) / p0 - 1
        gamma0 = E0 / mass
        beta0 = math.sqrt(1 - 1.0 / gamma0 / gamma0)

        if update_eta:
            delta_gamma = dEl / mass
        else:
            delta_gamma =0
        eta = alphac - 1.0/(gamma0+delta_gamma)/(gamma0+delta_gamma)


        pl += 2.0 * np.pi * harm * eta * dl
        phi_list.append(pl)
        dE_list.append(dEl)
        delta_list.append(dl)

    return np.vstack(phi_list), np.vstack(dE_list), np.vstack(delta_list)


'''

E00=100e9 
mass=938e6
gamma=E00/mass
beta=pylab.sqrt(1-1.0/gamma/gamma)

eV=5e6
alphac=0.002
gammat=pylab.sqrt(1/alphac)

h=360
ita=alphac-1/(E00*E00/mass/mass)
ita0=ita
print(ita)
nus=pylab.sqrt(h*abs(ita)*eV/2/pylab.pi/E00/beta/beta)
print(nus)


turns=5000
#inideltas=[0.001, 0.003, 0.008, 0.013, 0.017, 0.0175]
#iniphis=[phi_s, phi_s, phi_s, phi_s, phi_s, phi_s]
inideltas=[0.003]
iniphis=[phi_s, phi_s, phi_s,-pylab.pi*0.5,pylab.pi-phi_s]
#iniphis=[phi_s, phi_s, phi_s,pylab.pi*1.5,pylab.pi-phi_s]
for id in range(len(inideltas)):
    ini_delta=inideltas[id]
    ini_phi=iniphis[id]
    deltalist=[ini_delta,]
    philist=[ini_phi,]
    plotlist=[0,]
    E0=E00
    mass=938e6
    gamma=E0/mass
    beta=pylab.sqrt(1-1.0/gamma/gamma)
    beta0=beta
    eV=5e6
    relative_omega=1
    alphac=0.002
    gammat=pylab.sqrt(1/alphac)

    h=360
    ita=alphac-1/(E0*E0/mass/mass)
    for i in range(turns):
        dE0=deltalist[-1]*E0*beta*beta
        #dE0=deltalist[-1]
        phi0=philist[-1]
    
        dE1=dE0+eV*(pylab.sin(phi0)-pylab.sin(phi_s))
        oldE0=E0
        E0+=eV*pylab.sin(phi_s)
        
        delta_omega=-ita*relative_omega*(E0-oldE0)/oldE0/beta/beta
        relative_omega+=delta_omega
        gamma=E0/mass
        beta=pylab.sqrt(1-1.0/gamma/gamma)
        delta1=dE1/(E0*beta*beta)
        
        ita=alphac-1/gamma/gamma
        phi1=phi0+2*pylab.pi*h*ita*delta1
        plotlist.append(dE1/relative_omega)
        deltalist.append(delta1)
        philist.append(phi1)

    pylab.plot(philist[1:],plotlist[1:])
    pylab.plot(philist[0:1],deltalist[0:1],'r+')


pylab.xlabel("phase")
pylab.ylabel("energy deviation")
pylab.title('Turn {}'.format(turns))
#pylab.xlim([-pylab.pi/2,1.5*pylab.pi])
#pylab.ylim([-0.02,0.02])
pylab.show()
    
'''
