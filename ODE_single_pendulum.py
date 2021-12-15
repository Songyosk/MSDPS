"""
Module to solve ODE (i.e. equations of motion of a single pendulum) using various numerical methods

Author: Son Gyo Jung
Email: sgj30@cam.ac.uk
"""

from matplotlib import pyplot as plt
import numpy as np


class single_pendulum():
    """
    Class to solve equations of motion of a single pendulum using various numerical methods
    
    Note: Since the equations for the single pendulum motion are re-scaled into natural units, all the variables are in natural units.

    args:
        (1) h (type:float) - step size 
        (2) D (type:float) - damping coefficient
        (3) t - time steps defined using h

    """

    def __init__(self, h = 0.01, D = 0.):
        self.h = h
        self.D = D
        self.t = np.arange(0, 100, self.h)



    def EulerForward(self, save_figure = False):
        """
        Solve the ODE using Explicit Euler method

        args:
            (1) figure (type:bool) - whether to save the figure
            (2) w - angular velocity (starting from rest)
            (3) theta = angular displacement

        return:
            (1) Figure of angular displacement against time
        """    

        # Angular velocity (starting from rest)
        self.EulerForward_w = np.zeros(self.t.size) #starting from rest 
        
        # Angular displacement
        self.EulerForward_theta = np.ones(self.t.size) * 10. * np.pi / 180.  #* 0.75 * np.pi #for ExactEuler
        
        # Energy 
        self.EulerForward_energy = np.ones(self.t.size)
        

        for i in range(1, self.t.size):
            # Solve ODE using Euler Forward
            self.EulerForward_w[i] = self.EulerForward_w[i - 1] + (-self.D * self.EulerForward_w[i - 1] - self.EulerForward_theta[i - 1]) * self.h

            self.EulerForward_theta[i] = self.EulerForward_theta[i - 1] + self.EulerForward_w[i - 1] * self.h

            self.EulerForward_energy[i] = (0.5 * self.EulerForward_w[i]**(2)) + (0.5 * self.EulerForward_theta[i]**(2))
        

        plt.figure(figsize=(6, 4), dpi = 300)    

        plt.plot(self.t, self.EulerForward_theta, '.-', label = 'Explicit Euler with h=' + str(self.h) + ' , $\hat{D}$=' + str(self.D))
        plt.xlabel('t('+r'$\sqrt{\frac{l}{g}}$'+')', fontsize=14)
        plt.ylabel('$\Theta$ (rad)', fontsize=14)
        plt.legend(loc=1, prop={'size':14})
        plt.grid(True)
        plt.show()

        if save_figure == True:
            plt.savefig('Angular_displacement_EulerForward.png', dpi = 300, bbox_inches = "tight")




    def Leapfrog(self, save_figure = False):       
        """
        Solve the ODE using Leapfrog or Leapfrog with Euler Forward
        """
        EulerForwardPrediction = str.lower(input('Predict using Euler Forward for Leapfrog? "yes" or "no": '))    

        self.Leapfrog_w = np.zeros(self.t.size) #starting from rest 
        self.Leapfrog_theta = np.ones(self.t.size) * 10. * np.pi / 180.
                

        if EulerForwardPrediction == 'no':    
            for i in range(1, self.t.size):
                self.Leapfrog_w[i] = self.Leapfrog_w[i - 2] + 2 * (-self.D * self.Leapfrog_w[i - 1] - self.Leapfrog_theta[i - 1]) * self.h
                self.Leapfrog_theta[i] = self.Leapfrog_theta[i - 2] + 2 * self.Leapfrog_w[i - 1] * self.h 
            

            plt.figure(figsize=(6, 4), dpi = 300) 

            plt.plot(self.t, self.Leapfrog_theta, '.-', label = 'Leapfrog with h=' + str(self.h) + ' , $\hat{D}$=' + str(self.D))
            plt.xlabel('t('+r'$\sqrt{\frac{l}{g}}$'+')')
            plt.ylabel('$\Theta$ (rad)')
            plt.legend(loc=1, prop={'size':14})
            plt.title('Angular displacement vs time', fontsize=14)
            plt.grid(True)
            plt.show()


            if save_figure == True:
                plt.savefig('Angular_displacement_Leapfrog.png', dpi = 300, bbox_inches = "tight")
        

        elif EulerForwardPrediction == 'yes':   
            for i in range(2, self.t.size): #range 2 if predicting using EulerForward 
                self.Leapfrog_w[1] = self.Leapfrog_w[0] + (-self.D * self.Leapfrog_w[0] - self.Leapfrog_theta[0]) * self.h #Predict only the next value required for leapfrog method
                self.Leapfrog_theta[1] = self.Leapfrog_theta[0] + self.Leapfrog_w[0] * self.h

                self.Leapfrog_w[i] = self.Leapfrog_w[i - 2] + 2 * (-self.D * self.Leapfrog_w[i - 1] - self.Leapfrog_theta[i - 1]) * self.h
                self.Leapfrog_theta[i] = self.Leapfrog_theta[i - 2] + 2 * self.Leapfrog_w[i - 1] * self.h  
                

            plt.figure(figsize=(6, 4), dpi = 300)   

            plt.plot(self.t, self.Leapfrog_theta, '.-', label = 'Leapfrog with h=' + str(self.h) + ' , $\hat{D}$=' + str(self.D))
            plt.xlabel('t('+r'$\sqrt{\frac{l}{g}}$'+')')
            plt.ylabel('$\Theta$ (rad)')
            plt.legend(loc=1, prop={'size':14})
            plt.title('Angular displacement vs time', fontsize=14)
            plt.grid(True)
            plt.show()


            if save_figure == True:
                plt.savefig('Angular_displacement_Leapfrog_EulerForward.png', dpi = 300, bbox_inches = "tight")
        
            
        else:
            return('Invalid input')
            
        
        
 
    def EulerBackward(self, save_figure = False):
        """
        Solve the ODE using Implicit Euler
        """

        self.EulerBackward_w = np.zeros(self.t.size)
        self.EulerBackward_theta = np.ones(self.t.size) * 10. * np.pi / 180.
        

        for i in range(1, self.t.size):
            self.EulerBackward_w[i] = (-1 * self.h * self.EulerBackward_theta[i - 1] + self.EulerBackward_w[i - 1]) / (1 + self.h * self.D + self.h**(2))
            self.EulerBackward_theta[i] = ((1 + self.h * self.D) * self.EulerBackward_theta[i - 1] + self.EulerBackward_w[i] * self.h) / (1 + self.h * self.D + self.h**(2))    
        

        plt.figure(figsize=(6, 4), dpi = 300) 

        plt.plot(self.t, self.EulerBackward_theta, 'k.-', label = 'Implicit Euler with h=' + str(self.h) + ' , $\hat{D}$=' + str(self.D))
        plt.xlabel('t('+r'$\sqrt{\frac{l}{g}}$'+')')
        plt.ylabel('$\Theta$ (rad)')
        plt.legend(loc=1, prop={'size':14})
        plt.title('Angular displacement vs time', fontsize=22)
        plt.grid(True)
        plt.show()


        if save_figure == True:
            plt.savefig('Angular_displacement_EulerBackward.png', dpi = 300, bbox_inches = "tight")
        



    def RK4(self, save_figure = False):  
        """
        Solve the ODE using RK4
        """     

        self.RK4_w = np.zeros(self.t.size)
        self.RK4_theta = np.ones(self.t.size) * 10. * np.pi / 180.
        
        for i in range(1, self.t.size):
            f_w1 = -self.D * self.RK4_w[i - 1] - self.RK4_theta[i - 1]
            f_theta1 = self.RK4_w[i - 1]
            w1 = self.RK4_w[i - 1] + f_w1 * self.h/2.
            theta1 = self.RK4_theta[i - 1] + f_theta1 * self.h / 2.
            
            f_w2 = -self.D * w1 - theta1
            f_theta2 = w1
            w2 = self.RK4_w[i - 1] + f_w2 * self.h / 2.
            theta2 = self.RK4_theta[i - 1] + f_theta2 * self.h / 2.
            
            f_w3 = -self.D * w2 - theta2
            f_theta3 = w2
            w3 = self.RK4_w[i-1] + f_w3 * self.h
            theta3 = self.RK4_theta[i-1] + f_theta3 * self.h
            
            f_w4 = -self.D * w3 - theta3
            f_theta4 = w3
            
            self.RK4_w[i] = self.RK4_w[i - 1] + (f_w1 + 2 * f_w2 + 2 * f_w3 + f_w4)/6. * self.h
            self.RK4_theta[i] = self.RK4_theta[i - 1] + (f_theta1 + 2 * f_theta2 + 2 * f_theta3 + f_theta4) / 6 * self.h

        
        plt.figure(figsize=(6, 4), dpi = 300)   

        plt.plot(self.t, self.RK4_theta, 'r--', linewidth = 2, label = 'RK4 with h=' + str(self.h) + ' , $\hat{D}$=' + str(self.D))
        plt.xlabel('t('+r'$\sqrt{\frac{l}{g}}$'+')')
        plt.ylabel('$\Theta$ (rad)')
        plt.legend(loc=1, prop={'size':14})
        plt.title('Angular displacement vs time', fontsize=14)
        plt.grid(True)
        plt.show()


        if save_figure == True:
            plt.savefig('Angular_displacement_RK4.png', dpi = 300, bbox_inches = "tight")
        



    def Energy(self, save_figure = False, run_all = False):
        """
        Function to plot the total energy against time
        """
        
        if run_all == True:
            self.EulerForward()
            self.Leapfrog()
            self.EulerBackward()
            self.RK4()


        Energy1, Energy2, Energy3, Energy4 = [], [], [], []

        for i in range(self.t.size):
            Energy1.append(0.5 * self.EulerForward_w[i]**(2) + 0.5 * self.EulerForward_theta[i]**(2))
            Energy2.append(0.5 * self.Leapfrog_w[i]**(2) + 0.5 * self.Leapfrog_theta[i]**(2))
            Energy3.append(0.5 * self.EulerBackward_w[i]**(2) + 0.5 * self.EulerBackward_theta[i]**(2))
            Energy4.append(0.5 * self.RK4_w[i]**(2) + 0.5 * self.RK4_theta[i]**(2))
                 
                    
        plt.figure(figsize=(6, 4), dpi = 300)  

        plt.plot(self.t, Energy1, '-', linewidth = 2, label='Explicit Euler with h=' +str(self.h) + ' , $\hat{D}$=' + str(self.D))
        plt.plot(self.t, Energy2, '-', linewidth = 2, label='Leapfrog with h=' +str(self.h) + ' , $\hat{D}$=' + str(self.D))
        plt.plot(self.t, Energy3, 'k',linewidth = 2, label='Implicit Euler with h=' +str(self.h) + ' , $\hat{D}$=' + str(self.D))
        plt.plot(self.t, Energy4, 'r--',linewidth = 2, label='RK4 with h=' +str(self.h) + ' , $\hat{D}$=' + str(self.D))

        plt.legend(loc=1,prop={'size':14})
        plt.xlabel('t('+r'$\sqrt{\frac{l}{g}}$'+')', fontsize=18)
        plt.ylabel('E (mgl)', fontsize=18)
        plt.title('Total energy vs time', fontsize=18)
        plt.grid(True)
        plt.show()


        if save_figure == True:
            plt.savefig('Energy.png', dpi = 300, bbox_inches = "tight")
        



    def ExactEuler(self, save_figure = False, run_EulerForward = True, h = 0.3, D = 0.2): 
        """
        Function to plot the dynamics predicted with and without the small angle approximation

        Note: Don't forget to change the angle in EulerForward function as well as adjusting h and D
        """

        self.D = D
        self.h = h
        self.t = np.arange(0, 100, self.h)


        w = np.zeros(self.t.size) #initial speed = 0      
        theta = np.ones(self.t.size) * 0.75 * np.pi
        energy = np.zeros(self.t.size)


        if run_EulerForward == True:
            self.EulerForward()


        for i in range(1, self.t.size):
            w[i] = w[i - 1] + (-self.D * w[i-1] - np.sin(theta[i-1])) * self.h
            theta[i] = theta[i - 1] + w[i - 1] * self.h
            energy[i] = 0.5*(w[i]**2) + 1 - np.cos(theta[i])
        


        fig = plt.figure(figsize = (15, 10))

        ax1 = fig.add_subplot(211)
        ax1.plot(self.t, theta, label = 'Without small angle approx: h=' + str(self.h) + ' , $\hat{D}$=' + str(self.D))
        ax1.plot(self.t, self.EulerForward_theta, label = 'With small angle approx: h=' + str(self.h) + ' , $\hat{D}$=' + str(self.D))

        plt.legend() 
        plt.grid(True)
        plt.legend(loc=2,prop={'size':14})
        plt.xlabel('t('+r'$\sqrt{\frac{l}{g}}$'+')', fontsize=14)
        plt.ylabel('$\Theta$ (rad)', fontsize=14)
        ax1.set_title('Angular displacement vs time using explicit Euler for large initial angle', fontsize=14)
        plt.tight_layout()


        ax2 = fig.add_subplot(212)  
        ax2.plot(self.t, energy, label = 'Without small angle approx: h=' + str(self.h) + ' , $\hat{D}$=' + str(self.D))
        ax2.plot(self.t, self.EulerForward_energy, label = 'With small angle approx: with h=' + str(self.h) + ' , $\hat{D}$=' + str(self.D))

        plt.legend()
        plt.grid(True)
        plt.legend(loc=2,prop={'size':14})
        plt.xlabel('t('+r'$\sqrt{\frac{l}{g}}$'+')', fontsize=14)
        plt.ylabel('E (mgl)', fontsize=14)
        ax2.set_title('Total energy vs time using explicit Euler for large initial angle', fontsize=14)
        plt.tight_layout()
        plt.show()

        if save_figure == True:
            plt.savefig('compare.png', dpi = 300, bbox_inches = "tight")


        

    def period(self, method = 'EulerForward'):
        """
        Calculating the period of the oscillations
        """

        if method == 'EulerForward':
            method_theta = self.EulerForward_theta

        elif method == 'EulerBackward':
            method_theta = self.EulerBackward_theta

        elif method == 'Leapfrog':
            method_theta = self.Leapfrog_theta

        elif method == 'RK4':
            method_theta = self.RK4_theta


        period_xintercept = []
        period_time = []

        for i in range(1, self.t.size):
            if abs(method_theta[i]) <= 0.01:
                if method_theta[i] < 0 and method_theta[i + 1] > 0 or method_theta[i] > 0 and method_theta[i + 1] < 0:
                    period_xintercept.append(round(self.t[i], 3))
                    period_time.append(i)
        
        print('x-intercept at: \n', period_xintercept, '\n')  

        
        period_difference = []   

        for i in range(len(period_xintercept) - 1):
                period_difference.append(period_xintercept[i + 1] - period_xintercept[i])   

        print('period difference: \n', period_difference)            
                                        
