# Main Functions
## attitude control
Under Developement.

# Test Class
## Optimal-Control-Dynamic-Programming
MATLAB Optimal Control codes related to HJB Dynamic Programming to find the optimal path for any state of a linear system

## System with state equation A*X + B*U
This code solves the example at the end of chapter 3 of Optimal Control Theory (kirk)

## How to Run
The class is loaded with a example. In order to see the results, create a Dynamic_Solver object first:  
`objA = Dynamic_Solver;`  
Calculate the optimal controls at each state:  
`run(objA)`  
Finally, track the optimal states and controls:  
`get_optimal_path(objA)`

## Unique Performance
This code exploits Vectorization in MATLAB thus it is more than thousands orders of magnitude faster than loop algorithm implementation which is explained in Optimal Control Theory book(s).
