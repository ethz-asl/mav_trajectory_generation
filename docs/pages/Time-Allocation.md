# Time Allocation

Optimization for:
* Time only ([Burri](#burri-et-al), [Richter](#richter-et-al), [Mellinger](#melling-and-kumar) and [Segment violation](#segment-violation))
* Time and free constraints ([Burri](#burri-et-al) and [Richter](#richter-et-al))

Methods implemented:  
1. C. Richter, A. Bry, and N. Roy, “**Polynomial trajectory planning for aggressive quadrotor flight in dense indoor environments,**” in *International Journal of Robotics Research*, Springer, 2016.
2. M. Burri, H. Oleynikova, M. Achtelik, and R. Siegwart, “**Real-Time Visual-Inertial Mapping, Re-localization and Planning Onboard MAVs in Previously Unknown Environments**”. In *IEEE Int. Conf. on Intelligent Robots and Systems* (IROS), September 2015.
3. D. Mellinger and V. Kumar, "**Minimum Snap Trajectory Generation and Control for Quadrotors**"
4. Segment violation  


### Benchmark
* trajectory time
* computation time
* relative violation of velocity
* maximum distance between trajectory and straight line
* area between trajectory and straight line

Additionally:
* comparison of convergence time and quality of default and custom initial step
* comparison magic fabian vs. trapezoidal for initial time segments


### [Richter et al.](#richter-et-al)
Paper: **Polynomial trajectory planning for aggressive quadrotor flight in dense indoor environments**  
Published in: *International Journal of Robotics Research*, Springer  
Year: 2016   

Usable for optimization of: 
* Time only (`NonlinearOptimizationParameters::kRichterTime`):    
* Time and free derivatives (`NonlinearOptimizationParameters::kRichterTimeAndConstraints`):  
<img src="https://user-images.githubusercontent.com/17544220/38020829-8bba3e8e-327b-11e8-9b36-5ca2cb3fb609.png" height="100" />  <br />


### [Burri et al.](#burri-et-al)
Paper: **Real-Time Visual-Inertial Mapping, Re-localization and Planning Onboard MAVs in Previously Unknown Environments**  
Published in: *IEEE Int. Conf. on Intelligent Robots and Systems* (IROS)  
Year: 2015  

Usable for optimization of: 
* Time only (`NonlinearOptimizationParameters::kSquaredTime`)
* Time and free derivatives (`NonlinearOptimizationParameters::kSquaredTimeAndConstraints`):  
<img src="https://user-images.githubusercontent.com/17544220/38019910-2e6f41c2-3279-11e8-9602-1e861b05b520.png" height="100" />  <br />
<img src="https://user-images.githubusercontent.com/17544220/38019908-2e30f9bc-3279-11e8-8fdc-5f3d997da955.png" height="50" />  <br />
<img src="https://user-images.githubusercontent.com/17544220/38019909-2e538e3c-3279-11e8-9dee-fe08a881743c.png" height="80" />  <br />


### [Mellinger and Kumar](#mellinger-and-kumar)
Paper: **Minimum Snap Trajectory Generation and Control for Quadrotors**    
Published in: *IEEE International Conference on Robotics and Automation* (ICRA)  
Year: 2011  

Usable for optimization of: 
* Time only (`NonlinearOptimizationParameters::kMellingerOuterLoop`):  
<img src="https://user-images.githubusercontent.com/17544220/38021405-0544f040-327d-11e8-9d92-86d5c0c31f70.png" height="30" />  <br />
<img src="https://user-images.githubusercontent.com/17544220/38021406-0560c428-327d-11e8-87c6-8d0b1ccf1886.png" height="60" />  <br />
<img src="https://user-images.githubusercontent.com/17544220/38021407-0583bf46-327d-11e8-8ac5-d784736da1ad.png" height="100" />  <br />  


### [Segment violation](#segment-violation)

Usable for optimization of: 
* Time only:  
