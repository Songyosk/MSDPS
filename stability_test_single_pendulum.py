'''
Son-Gyo Jung 

The FDMs are identical to the Project A - Single.py.
This is to calculate and plot results corresponding to the critical step size.
'''

from matplotlib import pyplot as plt
import numpy as np

Global = str.lower(raw_input('Do you want D = 0? Yes or no: '))

while Global is not 'yes' or 'no':
    if Global == 'yes':
        D = 0.
        break 
    elif Global == 'no':
        D = float(raw_input('D:')) #Damping
        break
    else:
       Global = str.lower(raw_input('Do you want D = 0? Yes or no: '))


        
def EulerForward(h):
    h = h
    t = np.arange(0,50,h)
        
    EulerForward.w = np.zeros(t.size) #starting from rest 
    EulerForward.theta = np.ones(t.size)*1.*np.pi/180.
    
    for i in range(1, t.size):
        EulerForward.w[i] = EulerForward.w[i-1] + (-D*EulerForward.w[i-1] - EulerForward.theta[i-1])*h
        EulerForward.theta[i] = EulerForward.theta[i-1] + EulerForward.w[i-1]*h
    
def Leapfrog(h):
    h = h
    t = np.arange(0,50,h)
    
            
    Leapfrog.w = np.zeros(t.size) #starting from rest 
    Leapfrog.theta = np.ones(t.size)*1.*np.pi/180. 
              
    for i in range(2, t.size): #range 2 if predicting using EulerForward 
        Leapfrog.w[1] = Leapfrog.w[0] + (-D*Leapfrog.w[0] - Leapfrog.theta[0])*h #Predict only the next value required for leapfrog method
        Leapfrog.theta[1] = Leapfrog.theta[0] + Leapfrog.w[0]*h
        Leapfrog.w[i] = Leapfrog.w[i-2] + 2*(-D*Leapfrog.w[i-1] - Leapfrog.theta[i-1])*h
        Leapfrog.theta[i] = Leapfrog.theta[i-2] + 2*Leapfrog.w[i-1]*h  
        

def EulerBackward(h):
    h = h
    t = np.arange(0,50,h)
    
    EulerBackward.w =np.zeros(t.size)
    EulerBackward.theta = np.ones(t.size)*1.*np.pi/180.
    
    for i in range(1, t.size):
        EulerBackward.w[i] = (-1 * h * EulerBackward.theta[i-1] + EulerBackward.w[i-1])/(1 + h*D + h**(2))
        EulerBackward.theta[i] = ((1 + h * D) * EulerBackward.theta[i-1] + EulerBackward.w[i]*h)/(1 + h*D + h**(2))    

def RK4(h):
    h = h
    t = np.arange(0,50,h)
    RK4.w = np.zeros(t.size)
    RK4.theta = np.ones(t.size)*1.*np.pi/180.
    
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


        
def Critical(method):
    '''
    Calculates the critical step size for each method. 
    Must write the correct method as an input to this function.
    ''' 
    Critical.Energy = []
    small_h = np.arange(0.00001, 1, 0.001)
    Critical.h_values = []
    last_theta = []
    last_w = []
    unstable_h = []
    for i in small_h:
        method(i)
        final_E = round((0.5 * method.w[-1]**(2) + 0.5 * method.theta[-1]**(2)), 10)
        initial_E = round((0.5 * method.w[0]**(2) + 0.5 * method.theta[0]**(2)), 10)
        last_theta.append(method.theta[-1])
        last_w.append(method.w[-1])
                                          
        if 0.9*abs(initial_E) <= final_E <= 1.1*abs(initial_E):
            Critical.Energy.append(final_E)
            Critical.h_values.append(i)
        else:
            unstable_h.append(i)
                
    print('Out of range when h = ' + str(unstable_h[0]))       
         
    plt.figure(0)            
    plt.plot(Critical.h_values, Critical.Energy, label=method.__name__)
    plt.xlabel('h('+r'$\sqrt{\frac{l}{g}}$'+')', fontsize=18)
    plt.ylabel('E (mgl)', fontsize=18)
    plt.title('Total energy vs step size h ($\hat D$=' + str(D) +')', fontsize=18)
    plt.legend(loc=2,prop={'size':18})
    plt.grid(True)
    plt.show()
    
    AllEnergy = []
    for i in range(small_h.size):
        AllEnergy.append(0.5 * last_theta[i]**(2) + 0.5 * last_w[i]**(2)) 
    plt.figure(1)
    plt.plot(small_h, AllEnergy, label=method.__name__)
    plt.xlabel('h('+r'$\sqrt{\frac{l}{g}}$'+')', fontsize=18)
    plt.ylabel('E (mgl)', fontsize=18)
    plt.title('Total energy vs step size h ($\hat D$=' + str(D) +')', fontsize=18)
    plt.legend(loc=2,prop={'size':18})
    plt.grid(True)
    plt.show()
    
    