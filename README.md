# Modelling Single & Double Pendulum Systems (MSDPS)
Using finite difference methods for solving orrdinary differential equations



## Introduction
The dynamics of a single pendulum and a double pendulum were computationally examined to analyse different finite difference methods for solving ordinary differential equations (ODEs). The fourth-order Runge-Kutta (RK4) was found to be the most suitable method considering the stability and the accuracy it can provide when solving general oscillatory problems. Subsequently, RK4 was employed to model the dynamics of a double pendulum. This yielded the understanding of the motion and the variation of the energy over time from which the dependence of the system on the initial conditions and the limitations of RK4 could be inferred.



## Requirements
Python 2.x is required to run the scripts (except those with name beginning with 'ODE_').

Create an environment using conda as follows:
```bash
  conda create -n python2 python=2.x
```
Then activate the new environment by:
```bash
  conda activate python2
```

Python 3.x is required to run the scripts whose name begins with 'ODE_'.

Create an environment using conda as follows:
```bash
  conda create -n python3 python=3.x
```
Then activate the new environment by:
```bash
  conda activate python3
```



## Results

- *Single Pendulum*

![FIG1AmpvsT4noDamping](https://user-images.githubusercontent.com/56391325/146282678-98a5eeef-1c46-42fd-8ca4-726388c4ad86.png)

Figure 1: Dynamics of an undamped single pendulum using different FDMs. This illustrates the oscillations predicted by Leapfrog and RK4 are periodic with relatively constant amplitude. In contrast, the oscillations are becoming unbound over time for explicit Euler, whereas it appears to be converging for implicit Euler


![FIG2EvsT4noDamping](https://user-images.githubusercontent.com/56391325/146282888-b93a859f-6389-4cfd-b84f-c0d76a4977f5.png)

Figure 2: The change in the total energy over time for an undamped single pendulum using the same parameters as in FIG. 1. Without damping, the law of energy conservation states that the total energy should be constant. However, the total energy is increasing and decreasing for explicit and implicit Euler, respectively.


![FIG3AmpvsT4Damping](https://user-images.githubusercontent.com/56391325/146283126-8251c4d3-c305-4a08-afbe-7d260b53e21b.png)
Figure 3: Dynamics of a damped (D = 0.2) single pendulum with parameters. This illustrates that the oscillations predicted by RK4, explicit and implicit Euler are converging with time. The oscillations model by the Leapfrog method appears to converge initially. However, it becomes unbound at higher t. 

![FIG4EvsT4Damping](https://user-images.githubusercontent.com/56391325/146283367-7aa24896-3265-44b3-8b05-333eaaa6112c.png)
Figure 4: The change in the total energy over time for a damped single pendulum using the same parameters as in FIG. 3. With damping, the law of energy conservation states that the total energy should decrease. This is true for FDMs except Leapfrog, where the total energy is decreasing at the beginning but increases at higher t.
