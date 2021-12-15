'''
Son-Gyo Jung 

Single pendulum:
Since the equations for the single pendulum motion are re-scaled into natural units, all the variables are in natural units

h = step size
D = Damping coefficient
t = time
w = angular velocity (starting from rest)
theta = angular displacement
'''

from matplotlib import pyplot as plt
import numpy as np


Global = str.lower(raw_input('Use global values? Yes or no: '))

while Global is not 'yes' or 'no':
    if Global == 'yes':
        h = 0.01 #Step-size
        D = 0.  #Damping
        t = np.arange(0, 100, h)
        break
    elif Global == 'no':
        h = float(raw_input('h:'))
        D = float(raw_input('D:'))
        t = np.arange(0, 100, h)
        break
    else:
        Global = str.lower(raw_input('Use global values? Yes or no: '))

def EulerForward():
    '''
    Explicit Euler
    '''      
    EulerForward.w = np.zeros(t.size) #starting from rest 
    EulerForward.theta = np.ones(t.size)*10.*np.pi/180.  #* 0.75 * np.pi #for ExactEuler
    EulerForward.energy = np.ones(t.size)
    
    for i in range(1, t.size):
        EulerForward.w[i] = EulerForward.w[i-1] + (-D*EulerForward.w[i-1] - EulerForward.theta[i-1])*h
        EulerForward.theta[i] = EulerForward.theta[i-1] + EulerForward.w[i-1]*h
        EulerForward.energy[i] = (0.5 * EulerForward.w[i]**(2)) + (0.5 * EulerForward.theta[i]**(2))
    
    plt.figure(1)     
    plt.plot(t,EulerForward.theta, '.-', label='Explicit Euler with h=' +str(h)+ ' , $\hat D$=' + str(D))
    plt.title('Angular displacement vs time', fontsize=22)
    plt.xlabel('t('+r'$\sqrt{\frac{l}{g}}$'+')', fontsize=22)
    plt.ylabel('$\Theta$ (rad)', fontsize=22)
    plt.legend(loc=1, prop={'size':14})
    plt.grid(True)
    plt.show()


def Leapfrog():       
    EulerForwardPrediction = str.lower(raw_input('Predict using Euler Forward for Leapfrog? Yes or No: '))    
    
    Leapfrog.w = np.zeros(t.size) #starting from rest 
    Leapfrog.theta = np.ones(t.size)*10.*np.pi/180.
            
    if EulerForwardPrediction == 'no':    
        for i in range(1, t.size):
            Leapfrog.w[i] = Leapfrog.w[i-2] + 2*(-D*Leapfrog.w[i-1] - Leapfrog.theta[i-1])*h
            Leapfrog.theta[i] = Leapfrog.theta[i-2] + 2*Leapfrog.w[i-1]*h 
        
        plt.figure(1)     
        plt.plot(t,Leapfrog.theta, '.-', label='Leapfrog with h=' +str(h)+ ' , $\hat D$=' + str(D))
        plt.xlabel('t('+r'$\sqrt{\frac{l}{g}}$'+')')
        plt.ylabel('$\Theta$ (rad)')
        plt.legend(loc=1, prop={'size':14})
        plt.title('Angular displacement vs time', fontsize=22)
        plt.grid(True)
        plt.show()
    
    elif EulerForwardPrediction == 'yes':   
        for i in range(2, t.size): #range 2 if predicting using EulerForward 
            Leapfrog.w[1] = Leapfrog.w[0] + (-D*Leapfrog.w[0] - Leapfrog.theta[0])*h #Predict only the next value required for leapfrog method
            Leapfrog.theta[1] = Leapfrog.theta[0] + Leapfrog.w[0]*h
            Leapfrog.w[i] = Leapfrog.w[i-2] + 2*(-D*Leapfrog.w[i-1] - Leapfrog.theta[i-1])*h
            Leapfrog.theta[i] = Leapfrog.theta[i-2] + 2*Leapfrog.w[i-1]*h  
            
        plt.figure(1)     
        plt.plot(t,Leapfrog.theta, '.-', label='Leapfrog with h=' +str(h) + ' , $\hat D$=' + str(D))
        plt.xlabel('t('+r'$\sqrt{\frac{l}{g}}$'+')')
        plt.ylabel('$\Theta$ (rad)')
        plt.legend(loc=1, prop={'size':14})
        plt.title('Angular displacement vs time', fontsize=22)
        plt.grid(True)
        plt.show()
        
    else:
        return('Invalid input')
            
        
               
 
def EulerBackward():
    '''
    Implicit Euler
    '''        
    EulerBackward.w =np.zeros(t.size)
    EulerBackward.theta = np.ones(t.size)*10.*np.pi/180.
    
    for i in range(1, t.size):
        EulerBackward.w[i] = (-1 * h * EulerBackward.theta[i-1] + EulerBackward.w[i-1])/(1 + h*D + h**(2))
        EulerBackward.theta[i] = ((1 + h * D) * EulerBackward.theta[i-1] + EulerBackward.w[i]*h)/(1 + h*D + h**(2))    
    
    plt.figure(1) 
    plt.plot(t,EulerBackward.theta, 'k.-', label='Implicit Euler with h=' +str(h) + ' , $\hat D$=' + str(D))
    plt.xlabel('t('+r'$\sqrt{\frac{l}{g}}$'+')')
    plt.ylabel('$\Theta$ (rad)')
    plt.legend(loc=1, prop={'size':14})
    plt.title('Angular displacement vs time', fontsize=22)
    plt.grid(True)
    plt.show()



def RK4():        
    RK4.w =np.zeros(t.size)
    RK4.theta = np.ones(t.size)*10.*np.pi/180.
    
    for i in range(1, t.size):
        f_w1 = -D * RK4.w[i-1] - RK4.theta[i-1]
        f_theta1 = RK4.w[i-1]
        w1 = RK4.w[i-1] + f_w1 * h/2.
        theta1 = RK4.theta[i-1] + f_theta1 * h/2.
        
        f_w2 = -D * w1 - theta1
        f_theta2 = w1
        w2 = RK4.w[i-1] + f_w2 * h/2.
        theta2 = RK4.theta[i-1] + f_theta2 * h/2.
        
        f_w3 = -D * w2 - theta2
        f_theta3 = w2
        w3 = RK4.w[i-1] + f_w3 * h
        theta3 = RK4.theta[i-1] + f_theta3 * h
        
        f_w4 = -D * w3 - theta3
        f_theta4 = w3
        
        RK4.w[i] = RK4.w[i-1] + (f_w1 + 2 * f_w2 + 2 * f_w3 + f_w4)/6. * h
        RK4.theta[i] = RK4.theta[i-1] + (f_theta1 + 2 * f_theta2 + 2 * f_theta3 + f_theta4)/6 * h
    
    plt.figure(1)      
    plt.plot(t,RK4.theta, 'r--', linewidth=2, label='RK4 with h=' +str(h) + ' , $\hat D$=' + str(D))
    plt.xlabel('t('+r'$\sqrt{\frac{l}{g}}$'+')')
    plt.ylabel('$\Theta$ (rad)')
    plt.legend(loc=1, prop={'size':22})
    plt.title('Angular displacement vs time', fontsize=22)
    plt.grid(True)
    plt.show()


def Energy():
    '''
    Function to plot the total energy against time.
    '''
    Energy1, Energy2, Energy3, Energy4 = [], [], [], []
    EulerForward()
    Leapfrog()
    EulerBackward()
    RK4()
    for i in range(t.size):
        Energy1.append(0.5 * EulerForward.w[i]**(2) + 0.5 * EulerForward.theta[i]**(2))
        Energy2.append(0.5 * Leapfrog.w[i]**(2) + 0.5 * Leapfrog.theta[i]**(2))
        Energy3.append(0.5 * EulerBackward.w[i]**(2) + 0.5 * EulerBackward.theta[i]**(2))
        Energy4.append(0.5 * RK4.w[i]**(2) + 0.5 * RK4.theta[i]**(2))
             
                
    plt.figure(0) 
    plt.plot(t, Energy1, '-', linewidth=2, label='Explicit Euler with h=' +str(h) + ' , $\hat D$=' + str(D))
    plt.plot(t, Energy2, '-', linewidth=2, label='Leapfrog with h=' +str(h) + ' , $\hat D$=' + str(D))
    plt.plot(t, Energy3, 'k',linewidth=2, label='Implicit Euler with h=' +str(h) + ' , $\hat D$=' + str(D))
    plt.plot(t, Energy4, 'r--',linewidth=2, label='RK4 with h=' +str(h) + ' , $\hat D$=' + str(D))
    plt.legend(loc=1,prop={'size':14})
    plt.xlabel('t('+r'$\sqrt{\frac{l}{g}}$'+')', fontsize=18)
    plt.ylabel('E (mgl)', fontsize=18)
    plt.title('Total energy vs time', fontsize=18)
    #axes = plt.gca()
    #axes.set_ylim([0.005,0.025])
    #axes.set_xlim([0.,40.])
    #axes.set_yscale('log')
    #axes.set_xscale('log')
    plt.grid(True)
    plt.show()


def ExactEuler(): 
    '''
    Function to plot the dynamics predicted with and without the small angle approximation.
    Don't forget to change the angle in EulerForward function as well as adjusting h and D.
    '''
    w = np.zeros(t.size) #initial speed=0      
    theta = np.ones(t.size) * 0.75 *np.pi
    Energy = np.zeros(t.size)
    EulerForward()
    for i in range(1, t.size):
        w[i] = w[i-1] + (- D*w[i-1] - np.sin(theta[i-1]))*h
        theta[i] = theta[i-1] + w[i-1]*h
        Energy[i] = 0.5*(w[i]**2) + 1 - np.cos(theta[i])
    
    fig = plt.figure() 
    
    ax1 = fig.add_subplot(211)
    #plt.figure(2)    
    ax1.plot(t, theta, label='Without small angle approx: h=' + str(h) + ' , $\hat D$=' + str(D))
    ax1.plot(t, EulerForward.theta, label = 'With small angle approx: h=' + str(h) + ' , $\hat D$=' + str(D))
    plt.legend() 
    plt.grid(True)
    plt.legend(loc=2,prop={'size':18})
    plt.xlabel('t('+r'$\sqrt{\frac{l}{g}}$'+')', fontsize=18)
    plt.ylabel('$\Theta$ (rad)', fontsize=18)
    ax1.set_title('Angular displacement vs time using explicit Euler for large initial angle', fontsize=18)
    plt.tight_layout()
    plt.show() 
    
    ax2 = fig.add_subplot(212)
    #plt.figure(3)
    ax2.plot(t, Energy, label='Without small angle approx: h=' + str(h) + ' , $\hat D$=' + str(D))
    ax2.plot(t, EulerForward.energy, label='With small angle approx: with h=' + str(h) + ' , $\hat D$=' + str(D))
    plt.legend()
    plt.grid(True)
    plt.legend(loc=2,prop={'size':18})
    plt.xlabel('t('+r'$\sqrt{\frac{l}{g}}$'+')', fontsize=18)
    plt.ylabel('E (mgl)', fontsize=18)
    ax2.set_title('Total energy vs time using explicit Euler for large initial angle', fontsize=18)
    plt.tight_layout()
    plt.show()

    
def period(method):
    '''
    Calculating the period of the oscillations.
    '''
    method()
    period.xintercept = []
    period.time = []
    for i in range(1, t.size):
        if abs(method.theta[i]) <= 0.01:
            if method.theta[i] < 0 and method.theta[i+1] > 0 or  method.theta[i] > 0 and method.theta[i+1] < 0:
                period.xintercept.append(round(t[i],3))
                period.time.append(i)
    
    print(period.xintercept)    
    
    period.difference = []        
    for i in range(len(period.xintercept)):
        period.difference.append(period.xintercept[i] - period.xintercept[i-1])                
                                    
    plt.figure(0) 
    plt.plot(period.time, period.difference, '-', linewidth=2, label='Explicit Euler with h=' +str(h) + ' , $\hat D$= ' + str(D))
    plt.xlabel('time_hat')
    plt.ylabel('total_energy_hat')
    plt.legend(loc=2, prop={'size':12})
    plt.grid(True)
    plt.show()














