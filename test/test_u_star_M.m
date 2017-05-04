%% Initialization
D = Dynamic_Solver
run(D)

%% Visualize Optimum U(x1,x2) Matrix at every stage
D.plot_u_star()
%% find optimal path using one stage U matrix ( steady state U) 
D.get_optimal_path([2;1],'ssu',190) % bad
D.get_optimal_path([2;1],'ssu',180) % controllable but inefficient
D.get_optimal_path([2;1],'ssu',135) % near optimal control obtained at stage 1
D.get_optimal_path([2;1],'ssu',128) %
D.get_optimal_path([2;1],'ssu',127) %
D.get_optimal_path([2;1],'ssu',30) %

%% find optimal path by stepping through stage 1 to N
D.get_optimal_path 