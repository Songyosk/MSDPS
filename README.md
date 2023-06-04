# Modelling Single & Double Pendulum Systems (MSDPS)
Using finite difference methods for solving ordinary differential equations

[![DOI](https://zenodo.org/badge/438801270.svg)](https://zenodo.org/badge/latestdoi/438801270)


## Introduction
The dynamics of a single pendulum and a double pendulum were computationally examined to analyse different finite difference methods for solving ordinary differential equations (ODEs). The fourth-order Runge-Kutta (RK4) was found to be the most suitable method considering the stability and the accuracy it can provide when solving general oscillatory problems. Subsequently, RK4 was employed to model the dynamics of a double pendulum. This yielded the understanding of the motion and the variation of the energy over time from which the dependence of the system on the initial conditions and the limitations of RK4 could be inferred.



## Requirements
Python 2.x is required to run the scripts (except for those with name beginning with 'ODE_').

Create an environment using conda as follows:
```bash
  conda create -n python2 python=2.x
```
Then activate the new environment by:
```bash
  conda activate python2
```

Python 3.x is required to run the scripts whose name begins with 'ODE_'. These scripts are Python 3.x version of 'Project A - Single.py' and 'Project A - Double.py'.

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

 <br />

![FIG2EvsT4noDamping](https://user-images.githubusercontent.com/56391325/146282888-b93a859f-6389-4cfd-b84f-c0d76a4977f5.png)

Figure 2: The change in the total energy over time for an undamped single pendulum using the same parameters as in FIG. 1. Without damping, the law of energy conservation states that the total energy should be constant. However, the total energy is increasing and decreasing for explicit and implicit Euler, respectively.

 <br />

![FIG3AmpvsT4Damping](https://user-images.githubusercontent.com/56391325/146283126-8251c4d3-c305-4a08-afbe-7d260b53e21b.png)

Figure 3: Dynamics of a damped (D = 0.2) single pendulum with parameters. This illustrates that the oscillations predicted by RK4, explicit and implicit Euler are converging with time. The oscillations model by the Leapfrog method appears to converge initially. However, it becomes unbound at higher t. 

 <br />

![FIG4EvsT4Damping](https://user-images.githubusercontent.com/56391325/146283367-7aa24896-3265-44b3-8b05-333eaaa6112c.png)

Figure 4: The change in the total energy over time for a damped single pendulum using the same parameters as in FIG. 3. With damping, the law of energy conservation states that the total energy should decrease. This is true for FDMs except Leapfrog, where the total energy is decreasing at the beginning but increases at higher t.

 <br />

![FIG5ExactEuler](https://user-images.githubusercontent.com/56391325/146284060-310458cb-d7a0-4b0f-837d-0b91ed7dfd21.png)

Figure 5: The dynamics of the single pendulum and the variation of energy over time with and without the small angle approximation. The parameters are: h = 0.01, D = 0.2, theta_initial = (pi x 3)/4.

 <br />

![extra](https://user-images.githubusercontent.com/56391325/146284305-027faed5-22b7-4da4-982c-85ea0cec9861.png)

Figure 6: The dynamics of the single pendulum and the variation of energy over time with and without the small angle approximation. The parameters are: h = 0.3, D =0.2, theta_initial = (pi x 3)/4. The graph without the small angle approximation is stable, whereas the former is unstable.

 <br />
 
 
 
 - *Double Pendulum*

![FIG7DoubleMotionNoDamp](https://user-images.githubusercontent.com/56391325/146284589-d6e18c47-cdaa-45ac-811c-435abd73d7d0.png)

Figure 7: The dynamics of the double pendulum with parameters: h = 0.01, G = 0.0, R = {0.01, 1, 100}.

 <br />

![FIG9DoubleEnergyNoDamp](https://user-images.githubusercontent.com/56391325/146284688-47be8355-d5a7-497e-b05c-9ee539427597.png)

Figure 8: The variation of energy over time with parameters: h = 0.01, G = 0.0, R = {0.01, 1, 100} (i.e. the system in FIG. 7).

 <br />

![FIG8DoubleMotionDamp](https://user-images.githubusercontent.com/56391325/146284774-9bd5ab2a-1780-419b-8db9-290b99d115e3.png)

Figure 9: The dynamics of the double pendulum with parameters: h = 0.01, G = 1.0, R = {0.01, 1, 100}.

 <br />
 
![FIG10DoubleEnergyDamp](https://user-images.githubusercontent.com/56391325/146284872-c99fe9b6-7ae0-4f52-a1a0-717bc571d071.png)

Figure 10:  The variation of $E$ over $t$ with parameters: h = 0.01, G = 1.0, R = {0.01, 1, 100}.


## 🔗 Links
[![linkedin](https://img.shields.io/badge/S.G.Jung-0A66C2?style=for-the-badge&logo=linkedin&logoColor=white)](https://www.linkedin.com/in/son-gyo-jung-655537135/)


## License
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
