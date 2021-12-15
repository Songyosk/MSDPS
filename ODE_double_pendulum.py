"""
Module to solve ODE (i.e. equations of motion of a double pendulum) using various numerical methods

Author: Son Gyo Jung
Email: sgj30@cam.ac.uk
"""

import numpy as np
from matplotlib import pyplot as plt


class double_pendulum():
    """
    Class to solve equations of motion of a double pendulum 
    
    Note: Since the equations for the single pendulum motion are re-scaled into natural units, all the variables are in natural units.
    
    args:
        (1) h (type:float) - step size 
        (2) G (type:float) - a constant as defined in the paper
        (3) t - time steps defined using h
        (4) R = M/m (type:float) - a constant as defined in the paper
        (5) w & v (type:float) - angular velocities
        (6) theta & phi (type:float) - angular displacements

    """

    def __init__(self, h = 0.01, G = 0.):
        self.h = h
        self.G = G
        self.t = np.arange(0, 100, self.h)



    def extend(self, t_0, slope, h_value):
        """
        Function written to aid with the extrapolation method
        """
        t_1 = t_0 + slope * h_value

        return t_1



    def RK4(self, R):
        t = self.t
        G = self.G
        h = self.h

        RK4_theta = np.ones(t.size) * 0.1    
        RK4_w = np.zeros(t.size)           #d(theta)/dt
        RK4_phi = np.zeros(t.size) 
        RK4_v = np.zeros(t.size)       #d(phi)/dt
        RK4_Energy = np.zeros(t.size)
        

        for i in range(1, t.size):
            f_w1 = -1.*(R + 1.) * RK4_theta[i - 1] + R * RK4_phi[i - 1] - G * RK4_w[i - 1]
            f_theta1 = RK4_w[i - 1]
            f_v1 = (R + 1.) * RK4_theta[i - 1] - (R + 1.) * RK4_phi[i - 1] + G * (1. - 1./R) * RK4_w[i - 1] - G/R * RK4_v[i - 1]
            f_phi1 = RK4_v[i-1]
            
            f_w2 = -1. * (R + 1.) * self.extend(RK4_theta[i - 1], f_theta1, h/2.) + R * self.extend(RK4_phi[i - 1], f_phi1, h/2.) - G * self.extend(RK4_w[i - 1], f_w1, h/2.)
            f_theta2 = self.extend(RK4_w[i - 1], f_w1, h/2.)
            f_v2 = (R + 1.) * self.extend(RK4_theta[i - 1], f_theta1, h/2.) - (R + 1.) * self.extend(RK4_phi[i - 1], f_phi1, h/2.) + G * (1. - 1./R) * self.extend(RK4_w[i - 1], f_w1, h/2.) - G/R * self.extend(RK4_v[i - 1], f_v1, h/2.)
            f_phi2 = self.extend(RK4_v[i - 1], f_v1, h/2.)
            
            f_w3 = -1. * (R + 1.) * self.extend(RK4_theta[i - 1], f_theta2, h/2.) + R * self.extend(RK4_phi[i - 1], f_phi2, h/2.) - G * self.extend(RK4_w[i - 1], f_w2, h/2.)
            f_theta3 = self.extend(RK4_w[i - 1], f_w2, h/2.)
            f_v3 = (R + 1.) * self.extend(RK4_theta[i - 1], f_theta2, h/2.) - (R + 1.) * self.extend(RK4_phi[i - 1], f_phi2, h/2.) + G * (1. - 1./R) * self.extend(RK4_w[i - 1], f_w2, h/2.) - G/R * self.extend(RK4_v[i - 1], f_v2, h/2.)
            f_phi3 = self.extend(RK4_v[i - 1], f_v2, h/2.)
            
            f_w4 = -1. * (R + 1.) * self.extend(RK4_theta[i - 1], f_theta3, h) + R * self.extend(RK4_phi[i - 1], f_phi3, h) - G * self.extend(RK4_w[i - 1], f_w3, h)
            f_theta4 = self.extend(RK4_w[i - 1], f_w3, h)
            f_v4 = (R + 1.) * self.extend(RK4_theta[i - 1], f_theta3, h) - (R + 1.) * self.extend(RK4_phi[i - 1], f_phi3, h) + G * (1. - 1./R) * self.extend(RK4_w[i - 1], f_w3, h) - G/R * self.extend(RK4_v[i - 1], f_v3, h)
            f_phi4 = self.extend(RK4_v[i - 1], f_v3, h)
            
            
            RK4_w[i] = RK4_w[i - 1] + (f_w1 + 2 * f_w2 + 2 * f_w3 + f_w4 )/6. * h
            RK4_theta[i] = RK4_theta[i - 1] + (f_theta1 + 2 * f_theta2 + 2 * f_theta3 + f_theta4 )/6. * h
            RK4_v[i] = RK4_v[i - 1] + (f_v1 + 2 * f_v2 + 2 * f_v3 + f_v4 )/6. * h
            RK4_phi[i] = RK4_phi[i - 1] + (f_phi1 + 2 * f_phi2 + 2 * f_phi3 + f_phi4 )/6. * h


        return RK4_theta, RK4_w, RK4_phi, RK4_v, RK4_Energy



    def plotenergy(self):
        """
        Plots all the energy graphs for different parameters of R.
        """

        t = self.t
        G = self.G
        h = self.h

        Energy = [[], [], []]
        KE = [[], [], []]
        PE = [[], [], []] 
        R = [0.01, 1., 100.]


        for value in range(len(R)):
            RK4_theta, RK4_w, RK4_phi, RK4_v, RK4_Energy = self.RK4(R[value])
            #self.RK4(R[value])
            for i in range(t.size):
                KE[value].append(0.5 * ((R[value] + 1) * RK4_w[i]**(2) + R[value] * RK4_v[i]**(2) + 2 * R[value] * RK4_v[i] * RK4_w[i]))
                PE[value].append(0.5 * ((R[value] + 1) * RK4_theta[i]**(2) + R[value] * RK4_phi[i]**2))
                Energy[value].append(KE[value][i] + PE[value][i])


        fig = plt.figure(figsize = (15, 10)) 
        
        ax1 = fig.add_subplot(311)
        ax1.plot(t, KE[0], label='KE with h=' + str(h) + ' , G=' + str(G))
        ax1.plot(t, PE[0], label='PE with h=' + str(h) + ' , G=' + str(G))
        ax1.plot(t, Energy[0], label='Total E with h=' + str(h) + ' , G= ' + str(G))
        plt.xlabel('t('+r'$\sqrt{\frac{l}{g}}$'+')', fontsize=14)
        plt.ylabel('Energy (mgl)', fontsize=14)
        ax1.legend(loc=1,prop={'size':14})
        ax1.set_title('Energy vs time (R=0.01)', fontsize=14)
        plt.grid(True)
        plt.tight_layout()
        #plt.show()  
        
        #plt.figure(1)
        ax2 = fig.add_subplot(312)
        ax2.plot(t, KE[1], label='KE with h=' + str(h) + ' , G=' + str(G))
        ax2.plot(t, PE[1], label='PE with h=' + str(h) + ' , G=' + str(G))
        ax2.plot(t, Energy[1], label='Total E with h=' + str(h) + ' , G=' + str(G))
        plt.xlabel('t('+r'$\sqrt{\frac{l}{g}}$'+')', fontsize=14)
        plt.ylabel('Energy (mgl)', fontsize=14)
        ax2.legend(loc=1,prop={'size':14})
        ax2.set_title('Energy vs time (R=1)', fontsize=14)
        axes = plt.gca()
        axes.set_ylim([0.,0.011]) #Scale on y-axis
        plt.grid(True)
        plt.tight_layout()
        #plt.show()  
        
        #plt.figure(2)
        ax3 = fig.add_subplot(313)
        ax3.plot(t, KE[2], label='KE with h=' + str(h) + ' , G=' + str(G))
        ax3.plot(t, PE[2], label='PE with h=' + str(h) + ' , G=' + str(G))
        ax3.plot(t, Energy[2], label='Total E with h=' + str(h) + ' , G=' + str(G))
        plt.xlabel('t('+r'$\sqrt{\frac{l}{g}}$'+')', fontsize=14)
        plt.ylabel('Energy (mgl)', fontsize=14)
        ax3.set_title('Energy vs time (R=100)', fontsize=14)
        ax3.legend(loc=1,prop={'size':14})
        plt.grid(True)
        plt.tight_layout()
        plt.show()
                                                         
                                                                                                                                                                            
                                                                                 
    def plotmotion(self):
        '''
        Plots all the dynamics of a double pendulum for different values of R.
        '''

        t = self.t
        G = self.G
        h = self.h

        fig = plt.figure(figsize = (15, 10))  
          
        RK4_theta, RK4_w, RK4_phi, RK4_v, RK4_Energy = self.RK4(0.01)  
        ax1 = fig.add_subplot(311)
        ax1.plot(t, RK4_theta, label = 'theta $\Theta$ with h=' +str(h) + ' , G=' + str(G))
        ax1.plot(t, RK4_phi, label = 'phi $\Phi$ with h=' +str(h) + ' , G=' + str(G))
        plt.xlabel('t('+r'$\sqrt{\frac{l}{g}}$'+')', fontsize=18)
        plt.ylabel('Angular \n displacement (rad)', fontsize=15)
        

        RK4_theta, RK4_w, RK4_phi, RK4_v, RK4_Energy = self.RK4(1.)  
        ax2 = fig.add_subplot(312)
        ax2.plot(t, RK4_theta, label = 'theta $\Theta$ with h =' +str(h) + ' , G=' + str(G))
        ax2.plot(t, RK4_phi, label = 'phi $\Phi$ with h=' +str(h) + ' , G=' + str(G))
        plt.xlabel('t('+r'$\sqrt{\frac{l}{g}}$'+')', fontsize=18)
        plt.ylabel('Angular \n displacement (rad)', fontsize=15)
        
        RK4_theta, RK4_w, RK4_phi, RK4_v, RK4_Energy = self.RK4(100.)  
        ax3 = fig.add_subplot(313)
        ax3.plot(t, RK4_theta, label = 'theta $\Theta$ with h=' +str(h) + ' , G=' + str(G))
        ax3.plot(t, RK4_phi, label = 'phi $\Phi$ with h=' +str(h) + ' , G=' + str(G))
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
        
