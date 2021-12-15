'''
Son-Gyo Jung 

The FDMs are identical to the Project A - Double.py.
Bisection method to determine the critical step size which corresponds to the tolerance level epsilon.  
'''

from matplotlib import pyplot as plt
import numpy as np


D = 0.

def EulerForward(h, point):
    h = h

    EulerForward.t = np.arange(0,50,h)
        
    EulerForward.w = np.zeros(EulerForward.t.size) #starting from rest 
    EulerForward.theta = np.ones(EulerForward.t.size)*1.*np.pi/180.
    EulerForward.Energy = np.zeros(EulerForward.t.size)
    

    for i in range(1, EulerForward.t.size):
        EulerForward.w[i] = EulerForward.w[i-1] + (-D*EulerForward.w[i-1] - EulerForward.theta[i-1])*h
        EulerForward.theta[i] = EulerForward.theta[i-1] + EulerForward.w[i-1]*h  
    

    for i in range(EulerForward.t.size):    
        EulerForward.Energy[i] = (0.5 * EulerForward.w[i]**(2) + 0.5 * EulerForward.theta[i]**(2))
    

    if point == 'final':        
        return(round((0.5 * EulerForward.w[-1]**(2) + 0.5 * EulerForward.theta[-1]**(2)), 10))
    elif point == 'initial':
        return(round((0.5 * EulerForward.w[0]**(2) + 0.5 * EulerForward.theta[0]**(2)), 10))

        

def Leapfrog(h, point):
    h = h
    Leapfrog.t = np.arange(0,50,h) 
            
    Leapfrog.w = np.zeros(Leapfrog.t.size) #starting from rest 
    Leapfrog.theta = np.ones(Leapfrog.t.size)*1.*np.pi/180.
    Leapfrog.Energy = np.zeros(Leapfrog.t.size) 
              

    for i in range(2, Leapfrog.t.size): #range 2 if predicting using EulerForward 
        Leapfrog.w[1] = Leapfrog.w[0] + (-D*Leapfrog.w[0] - Leapfrog.theta[0])*h #Predict only the next value required for leapfrog method
        Leapfrog.theta[1] = Leapfrog.theta[0] + Leapfrog.w[0]*h
        Leapfrog.w[i] = Leapfrog.w[i-2] + 2*(-D*Leapfrog.w[i-1] - Leapfrog.theta[i-1])*h
        Leapfrog.theta[i] = Leapfrog.theta[i-2] + 2*Leapfrog.w[i-1]*h  
    

    for i in range(Leapfrog.t.size):   
        Leapfrog.Energy[i] = (0.5 * Leapfrog.w[i]**(2) + 0.5 * Leapfrog.theta[i]**(2))
        

    if point == 'final':     
        return(round((0.5 * Leapfrog.w[-1]**(2) + 0.5 * Leapfrog.theta[-1]**(2)), 10))

    elif point == 'initial':
        return(round((0.5 * Leapfrog.w[0]**(2) + 0.5 * Leapfrog.theta[0]**(2)), 10))

        
    
def EulerBackward(h, point):
    h = h
    EulerBackward.t = np.arange(0,50,h)
    
    EulerBackward.w =np.zeros(EulerBackward.t.size)
    EulerBackward.theta = np.ones(EulerBackward.t.size)*1.*np.pi/180.
    EulerBackward.Energy = np.zeros(EulerBackward.t.size) 
    

    for i in range(1, EulerBackward.t.size):
        EulerBackward.w[i] = (-1 * h * EulerBackward.theta[i-1] + EulerBackward.w[i-1])/(1 + h*D + h**(2))
        EulerBackward.theta[i] = ((1 + h * D) * EulerBackward.theta[i-1] + EulerBackward.w[i]*h)/(1 + h*D + h**(2))
    

    for i in range(EulerBackward.t.size):    
        EulerBackward.Energy[i] = (0.5 * EulerBackward.w[i]**(2) + 0.5 * EulerBackward.theta[i]**(2))
    

    if point == 'final':     
        return(round((0.5 * EulerBackward.w[-1]**(2) + 0.5 * EulerBackward.theta[-1]**(2)), 10)) 
    elif point == 'initial':
        return(round((0.5 * EulerBackward.w[0]**(2) + 0.5 * EulerBackward.theta[0]**(2)), 10)) 
           


def RK4(h, point):
    h = h
    RK4.t = np.arange(0,50,h)
    RK4.w = np.zeros(RK4.t.size)
    RK4.theta = np.ones(RK4.t.size)*1.*np.pi/180.
    RK4.Energy = np.zeros(RK4.t.size) 

    
    for i in range(1, RK4.t.size):
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
    

    for i in range(RK4.t.size):    
        RK4.Energy[i] = (0.5 * RK4.w[i]**(2) + 0.5 * RK4.theta[i]**(2))    
        

    if point == 'final':    
        return(round((0.5 * RK4.w[-1]**(2) + 0.5 * RK4.theta[-1]**(2)), 10))
        
    elif point == 'initial':
        return(round((0.5 * RK4.w[0]**(2) + 0.5 * RK4.theta[0]**(2)), 10))

        


def bisectionsearch():
    '''
    Bisection search method to estimate the critical step size.
    Must follow the instructions and submit a valid input. 
    Plots all the necessary results.
    '''
    method = None    
    input = str.lower(raw_input('Which method would you like to analyse: Eulerforward, Eulerbackward, leapfrog or RK4? '))
    while input is not 'eulerforward' or 'leapfrog' or 'eulerbackward' or 'rk4':
        if input == 'eulerforward':
            method = EulerForward          
            break
        elif input == 'leapfrog':
            method = Leapfrog
            break 
        elif input == 'eulerbackward':
            method = EulerBackward
            break 
        elif input == 'rk4':
            method = RK4
            break
        else:
            print('Invalid input! Please try again.')
            method = str.lower(raw_input('Which method would you like to analyse: Eulerforward, Eulerbackward, leapfrog or RK4? '))  
     
    h = float(raw_input('Take a guess for h: '))
    #epsilon = 0.000001
    numGuesses = 0
    low = 0.0001
    high = h
    avg = (high + low)/2.0
    epsilon = method(avg, 'initial')*0.1
    
    while abs(round(method(avg, 'final'), 10) - round(method(avg, 'initial'), 10)) >= epsilon:
        print('low h = ' + str(low) + ', high h = ' + str(high) + ', avg h = ' + str(avg))
        print('The energy difference is = ' +str(abs(method(avg, 'final') - method(avg, 'initial'))))
        numGuesses += 1
        
        if high == avg:
            print('Unstable FDM.')
            break
        elif low == avg:
            print('Unstable FDM.')
            break
            
        if method == EulerBackward:
            if method(avg, 'final') > method(avg, 'initial'):  
                low = avg
            else:
                high = avg
        else: 
            if method(avg, 'final') < method(avg, 'initial'):  
                low = avg
            else:
                high = avg
                 
        avg = (high + low)/2.0
     

    print('numGuesses = ' + str(numGuesses))
    print('h = ' + str(avg) + ' gives final energy = ' + str(method(avg, 'final')) + ' and initial energy = ' + str(method(avg, 'initial')))
    print('The energy difference is = ' +str(abs(method(avg, 'final') - method(avg, 'initial'))) + 'and the tolerance level is = ' + str(epsilon))
    
    
    
    fig = plt.figure() 
    #plt.figure(0)
    ax1 = fig.add_subplot(121)
      
    method(avg, 'a')
    ax1.plot(method.t, method.theta, label = 'theta $\Theta$ with h=' +str(avg) + ' , $\hat D$=' + str(D))
    plt.legend(loc=1, prop={'size':18})
    plt.xlabel('t('+r'$\sqrt{\frac{l}{g}}$'+')', fontsize=22)
    plt.ylabel('Angular \n displacement (rad)', fontsize=22)
    plt.grid(True)

    method(h, 'a')
    ax1.plot(method.t, method.theta, label = 'theta $\Theta$ with h=' +str(h) + ' , $\hat D$=' + str(D))
    plt.legend(loc=1, prop={'size':18})
    plt.xlabel('t('+r'$\sqrt{\frac{l}{g}}$'+')', fontsize=22)
    plt.ylabel('Angular \n displacement (rad)', fontsize=22)
    plt.grid(True)
    ax1.set_title('Angular displacement vs time (tolerance level: 10% of initial E)', fontsize=22)

    
    
    #plt.figure(1)
    ax2 = fig.add_subplot(122)
    
    method(avg, 'a')
    ax2.plot(method.t, method.Energy, label = 'Energy with h=' +str(avg) + ' , $\hat D$=' + str(D))
    plt.legend(loc=1, prop={'size':18})
    plt.xlabel('t('+r'$\sqrt{\frac{l}{g}}$'+')', fontsize=22)
    plt.ylabel('E (mgl)', fontsize=22)
    plt.grid(True)
    
    method(h, 'a')
    ax2.plot(method.t, method.Energy, label = 'Energy with h=' +str(h) + ' , $\hat D$=' + str(D))
    plt.legend(loc=1, prop={'size':18})
    plt.xlabel('t('+r'$\sqrt{\frac{l}{g}}$'+')', fontsize=22)
    plt.ylabel('E (mgl)', fontsize=22)
    plt.grid(True)
    ax2.set_title('Total energy vs time (tolerance level: 10% of initial E)', fontsize=22)
    plt.show()
