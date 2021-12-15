'''
Son-Gyo Jung 

Bisection method to determine the critical step size which corresponds to the tolerance level epsilon.
This was not used for this analysis but was created out of interest. 
Feel free to play around with it but can't guarantee it will yield any useful results.
'''     

from matplotlib import pyplot as plt
import numpy as np


   
G = 0.  #Damping

def extend(t_0, slope, h_value):
    t_1 = t_0 + slope * h_value
    return t_1



def RK4(R, h, point):
    RK4.t = np.arange(0,60,h)
    RK4.theta = np.ones(RK4.t.size) *0.1    
    RK4.w = np.zeros(RK4.t.size)           #d(theta)/dt
    RK4.phi = np.zeros(RK4.t.size) 
    RK4.v = np.zeros(RK4.t.size)       #d(phi)/dt
    RK4.Energy = np.zeros(RK4.t.size)
    
    for i in range(1, RK4.t.size):
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
    
    for i in range(RK4.t.size):
        RK4.Energy[i] = (0.5 * ((R+1)*RK4.w[i]**(2) + R*RK4.v[i]**(2) + 2*R*RK4.v[i]*RK4.w[i])) + (0.5 * ((R+1)*RK4.theta[i]**(2) + R*RK4.phi[i]**2))


    if point == 'final':
        return((0.5 * ((R+1)*RK4.w[-1]**(2) + R*RK4.v[-1]**(2) + 2*R*RK4.v[-1]*RK4.w[-1])) + (0.5 * ((R+1)*RK4.theta[-1]**(2) + R*RK4.phi[-1]**2)))

    elif point == 'initial':
        return((0.5 * ((R+1)*RK4.w[0]**(2) + R*RK4.v[0]**(2) + 2*R*RK4.v[0]*RK4.w[0])) + (0.5 * ((R+1)*RK4.theta[0]**(2) + R*RK4.phi[0]**2)))



def bisectionsearch():
       
    R = 0.01
    h = float(raw_input('Take a guess for h: '))
    epsilon = 0.00001
    numGuesses = 0
    low = 0.0001
    high = h
    avg = (high + low)/2.0
    
    while abs(round(RK4(R, avg, 'final'), 10) - round(RK4(R, avg, 'initial'), 10)) >= epsilon:
        print('low h = ' + str(low) + ', high h = ' + str(high) + ', avg h = ' + str(avg))
        numGuesses += 1
        
        if high == avg:
            print('Unstable FDM.')
            break
        elif low == avg:
            print('Unstable FDM.')
            break

         
        if RK4(R, avg, 'final') < RK4(R, avg, 'initial'):  
            low = avg
        else:
            high = avg
                 
        avg = (high + low)/2.0
     

    print('numGuesses = ' + str(numGuesses))
    print('h = ' + str(avg) + ' gives final energy = ' + str(RK4(R, avg, 'final')) + ' and initial energy = ' + str(RK4(R, avg, 'initial')))
    print('The energy difference is = ' +str(abs(RK4(R, avg, 'final') - RK4(R, avg, 'initial'))) + 'and the tolerance level is = ' + str(epsilon))
    
    fig = plt.figure() 
    #plt.figure(0) 
    ax1 = fig.add_subplot(121)
      
    RK4(R, avg, 'a')
    ax1.plot(RK4.t, RK4.theta, label = 'theta $\Theta$ with h =' +str(avg) + ' , D = ' + str(G))
    ax1.plot(RK4.t, RK4.phi, label = 'theta $\Theta$ with h =' +str(avg) + ' , D = ' + str(G))
    plt.legend(loc=1, prop={'size':10})
    plt.xlabel('time_hat')
    plt.ylabel('Angular \n displacement')
    plt.grid(True)

    RK4(R, h, 'a')
    ax1.plot(RK4.t, RK4.theta, label = 'theta $\Theta$ with h =' +str(h) + ' , D = ' + str(G))
    ax1.plot(RK4.t, RK4.phi, label = 'theta $\Theta$ with h =' +str(h) + ' , D = ' + str(G))
    plt.legend(loc=1, prop={'size':10})
    plt.xlabel('time_hat')
    plt.ylabel('Angular \n displacement')
    plt.grid(True)
    ax1.set_title('Angular displacement vs time (tolerance level: ' + str(epsilon) +')')

    
    
    #plt.figure(1)
    ax2 = fig.add_subplot(122)
    
    RK4(R, avg, 'a')
    ax2.plot(RK4.t, RK4.Energy, label = 'Energy with h =' +str(avg) + ' , D = ' + str(G))
    plt.legend(loc=1, prop={'size':10})
    plt.xlabel('time_hat')
    plt.ylabel('Total Energy')
    plt.grid(True)
    
    RK4(R, h, 'a')
    ax2.plot(RK4.t, RK4.Energy, label = 'Energy with h =' +str(4) + ' , D = ' + str(G))
    plt.legend(loc=1, prop={'size':10})
    plt.xlabel('time_hat')
    plt.ylabel('Total Energy')
    ax2.set_title('Energy vs time (tolerance level: ' + str(epsilon) +')')
    plt.grid(True)
    plt.show()
    
    
    