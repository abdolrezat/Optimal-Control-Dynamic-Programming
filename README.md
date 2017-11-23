
# Legacy Repository

This code repository contains codes for Dynamic Programming with application in Optimal Control for seperate position control and attitude control. It is now re-designed to incorporate full coupled position and attitude control of a satellite, see the successor repository: https://github.com/abdolrezat/SPHERES-DPcontrol 


## Run Test Class

MATLAB Optimal Control codes related to HJB Dynamic Programming to find the optimal path for any state of a linear system

The Test Class solves the example at the end of chapter 3 of Optimal Control Theory - kirk (System with state equation A*X + B*U
)

### Run the Tests
The class is loaded with a example. In order to see the results, create a Dynamic_Solver object first:  
`objA = Dynamic_Solver;`  
Calculate the optimal controls at each state:  
`run(objA)`  
Finally, track the optimal states and controls:  
`get_optimal_path(objA)`  
or get the optimal states and controls from any desired initial states:  
`get_optimal_path(objA, [x1 x2])`

### Unique Performance
This code exploits Vectorization in MATLAB thus it is more than thousands orders of magnitude faster than loop algorithm implementation which is explained in Optimal Control Theory book(s).
