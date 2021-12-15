'''
Son-Gyo Jung 

Double Pendulum
Since the equations for the double pendulum motion are re-scaled into natural units, all the variables are in natural units.

Global value by default

h = step size
G = a constant as defined in the paper
t = time
R = M/m #ratio between the masses
w, v = angular velocities 
theta, phi = angular displacements
'''


import numpy as np
from matplotlib import pyplot as plt


Global = str.lower(raw_input('Use global values? Yes or no: '))

while Global is not 'yes' or 'no':
    if Global == 'yes':
        h = 0.01 #Step-size
        G = 0.  #Damping
        t = np.arange(0, 100, h)
        break
    elif Global == 'no':
        h = float(raw_input('h:'))
        G = float(raw_input('G:'))
        t = np.arange(0, 100, h)
        break
    else:
        Global = str.lower(raw_input('Use global values? Yes or no: '))



def extend(t_0, slope, h_value):
    '''
    Function written to aid with the extrapolation method
    '''
    t_1 = t_0 + slope * h_value
    return t_1



def RK4(R):
    RK4.theta = np.ones(t.size) *0.1    
    RK4.w = np.zeros(t.size)           #d(theta)/dt
    RK4.phi = np.zeros(t.size) 
    RK4.v = np.zeros(t.size)       #d(phi)/dt
    RK4.Energy = np.zeros(t.size)
    
    for i in range(1, t.size):
        f_w1 = -1.*(R+1.)*RK4.theta[i-1] + R*RK4.phi[i-1] - G*RK4.w[i-1]
        f_theta1 = RK4.w[i-1]
        f_v1 = (R+1.)*RK4.theta[i-1] - (R+1.)*RK4.phi[i-1] + G*(1.-1./R)*RK4.w[i-1] - G/R * RK4.v[i-1]
        f_phi1 = RK4.v[i-1]
        
        f_w2 = -1.*(R+1.)*extend(RK4.theta[i-1], f_theta1, h/2.) + R*extend(RK4.phi[i-1], f_phi1, h/2.) - G*extend(RK4.w[i-1], f_w1, h/2.)
        f_theta2 = extend(RK4.w[i-1], f_w1, h/2.)
        f_v2 = (R+1.)*extend(RK4.theta[i-1], f_theta1, h/2.) - (R+1.)*extend(RK4.phi[i-1], f_phi1, h/2.) + G*(1.-1./R)*extend(RK4.w[i-1], f_w1, h/2.) - G/R * extend(RK4.v[i-1], f_v1, h/2.)
        f_phi2 = extend(RK4.v[i-1], f_v1, h/2.)
        
        f_w3 = -1.*(R+1.)*extend(RK4.theta[i-1], f_theta2, h/2.) + R*extend(RK4.phi[i-1], f_phi2, h/2.) - G*extend(RK4.w[i-1], f_w2, h/2.)
        f_theta3 = extend(RK4.w[i-1], f_w2, h/2.)
        f_v3 = (R+1.)*extend(RK4.theta[i-1], f_theta2, h/2.) - (R+1.)*extend(RK4.phi[i-1], f_phi2, h/2.) + G*(1.-1./R)*extend(RK4.w[i-1], f_w2, h/2.) - G/R * extend(RK4.v[i-1], f_v2, h/2.)
        f_phi3 = extend(RK4.v[i-1], f_v2, h/2.)
        
        f_w4 = -1.*(R+1.)*extend(RK4.theta[i-1], f_theta3, h) + R*extend(RK4.phi[i-1], f_phi3, h) - G*extend(RK4.w[i-1], f_w3, h)
        f_theta4 = extend(RK4.w[i-1], f_w3, h)
        f_v4 = (R+1.)*extend(RK4.theta[i-1], f_theta3, h) - (R+1.)*extend(RK4.phi[i-1], f_phi3, h) + G*(1.-1./R)*extend(RK4.w[i-1], f_w3, h) - G/R * extend(RK4.v[i-1], f_v3, h)
        f_phi4 = extend(RK4.v[i-1], f_v3, h)
        
        
        RK4.w[i] = RK4.w[i-1] + (f_w1 + 2 * f_w2 + 2 * f_w3 + f_w4 )/6. * h
        RK4.theta[i] = RK4.theta[i-1] + (f_theta1 + 2 * f_theta2 + 2 * f_theta3 + f_theta4 )/6. * h
        RK4.v[i] = RK4.v[i-1] + (f_v1 + 2 * f_v2 + 2 * f_v3 + f_v4 )/6. * h
        RK4.phi[i] = RK4.phi[i-1] + (f_phi1 + 2 * f_phi2 + 2 * f_phi3 + f_phi4 )/6. * h
        


def plotenergy():
    '''
    Plots all the energy graphs for different parameters of R.
    '''           
    Energy = [[], [], []]
    KE = [[], [], []]
    PE = [[], [], []] 
    R = [0.01, 1., 100.]

    for value in range(len(R)):
        RK4(R[value])
        for i in range(t.size):
            KE[value].append(0.5 * ((R[value]+1)*RK4.w[i]**(2) + R[value]*RK4.v[i]**(2) + 2*R[value]*RK4.v[i]*RK4.w[i]))
            PE[value].append(0.5 * ((R[value]+1)*RK4.theta[i]**(2) + R[value]*RK4.phi[i]**2))
            Energy[value].append(KE[value][i] + PE[value][i])
                    
    fig = plt.figure(0) 
    
    ax1 = fig.add_subplot(311)
    ax1.plot(t, KE[0], label='KE with h=' +str(h) + ' , G=' + str(G))
    ax1.plot(t, PE[0], label='PE with h=' +str(h) + ' , G=' + str(G))
    ax1.plot(t, Energy[0], label='Total E with h=' +str(h) + ' , G= ' + str(G))
    plt.xlabel('t('+r'$\sqrt{\frac{l}{g}}$'+')', fontsize=18)
    plt.ylabel('Energy (mgl)', fontsize=18)
    ax1.legend(loc=1,prop={'size':18})
    ax1.set_title('Energy vs time (R=0.01)', fontsize=18)
    plt.grid(True)
    plt.tight_layout()
    plt.show()  
    
    #plt.figure(1)
    ax2 = fig.add_subplot(312)
    ax2.plot(t, KE[1], label='KE with h=' +str(h) + ' , G=' + str(G))
    ax2.plot(t, PE[1], label='PE with h=' +str(h) + ' , G=' + str(G))
    ax2.plot(t, Energy[1], label='Total E with h=' +str(h) + ' , G=' + str(G))
    plt.xlabel('t('+r'$\sqrt{\frac{l}{g}}$'+')', fontsize=18)
    plt.ylabel('Energy (mgl)', fontsize=18)
    ax2.legend(loc=1,prop={'size':18})
    ax2.set_title('Energy vs time (R=1)', fontsize=18)
    axes = plt.gca()
    axes.set_ylim([0.,0.011]) #Scale on y-axis
    plt.grid(True)
    plt.tight_layout()
    plt.show()  
    
    #plt.figure(2)
    ax3 = fig.add_subplot(313)
    ax3.plot(t, KE[2], label='KE with h=' +str(h) + ' , G=' + str(G))
    ax3.plot(t, PE[2], label='PE with h=' +str(h) + ' , G=' + str(G))
    ax3.plot(t, Energy[2], label='Total E with h=' +str(h) + ' , G=' + str(G))
    plt.xlabel('t('+r'$\sqrt{\frac{l}{g}}$'+')', fontsize=18)
    plt.ylabel('Energy (mgl)', fontsize=18)
    ax3.set_title('Energy vs time (R=100)', fontsize=18)
    ax3.legend(loc=1,prop={'size':18})
    plt.grid(True)
    plt.tight_layout()
    plt.show()
                                                     
                                                                                                                                                                        
                                                                             
def plotmotion():
    '''
    Plots all the dynamics of a double pendulum for different values of R.
    '''
    fig = plt.figure()  
      
    RK4(0.01)
    ax1 = fig.add_subplot(311)
    ax1.plot(t, RK4.theta, label = 'theta $\Theta$ with h=' +str(h) + ' , G=' + str(G))
    ax1.plot(t, RK4.phi, label = 'phi $\Phi$ with h=' +str(h) + ' , G=' + str(G))
    plt.xlabel('t('+r'$\sqrt{\frac{l}{g}}$'+')', fontsize=18)
    plt.ylabel('Angular \n displacement (rad)', fontsize=15)
    

    RK4(1.)
    ax2 = fig.add_subplot(312)
    ax2.plot(t, RK4.theta, label = 'theta $\Theta$ with h =' +str(h) + ' , G=' + str(G))
    ax2.plot(t, RK4.phi, label = 'phi $\Phi$ with h=' +str(h) + ' , G=' + str(G))
    plt.xlabel('t('+r'$\sqrt{\frac{l}{g}}$'+')', fontsize=18)
    plt.ylabel('Angular \n displacement (rad)', fontsize=15)
    
    RK4(100.)
    ax3 = fig.add_subplot(313)
    ax3.plot(t, RK4.theta, label = 'theta $\Theta$ with h=' +str(h) + ' , G=' + str(G))
    ax3.plot(t, RK4.phi, label = 'phi $\Phi$ with h=' +str(h) + ' , G=' + str(G))
    plt.xlabel('t('+r'$\sqrt{\frac{l}{g}}$'+')', fontsize=18)
    plt.ylabel('Angular \n displacement (rad)', fontsize=15)

    plt.tight_layout()
    ax1.set_title('Angular displacement vs time (R=0.01)', fontsize=18)
    ax1.legend(loc=1,prop={'size':18})
    ax1.grid(True)
    
    ax2.set_title('Angular displacement vs time (R=1)', fontsize=18)
    ax2.legend(loc=1,prop={'size':18})
    ax2.grid(True)
    ax3.set_title('Angular displacement vs time (R=100)', fontsize=18)
    ax3.legend(loc=1,prop={'size':18})
    ax3.grid(True)
    
