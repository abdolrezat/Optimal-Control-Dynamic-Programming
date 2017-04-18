classdef Solver_attitude < dynamicprops
    %SOLVER_ATTITUDE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        N % number of stages
        J1
        J2
        J3
        dim_U1 % dimension number of U, after the 7 state variables
        dim_U2 
        dim_U3 
        U_vector % values that are applied as control inputs
        % big data
        J_next_states_opt single  % optimum values of J will be stored in this matrix and recalled each stage
        J_current_state_fix single  % values of cost to reach next stage, always fixed in Time Invariant system
        X1 single % (s+c)-dimensional grid of state 1
        X2 single 
        X3 single 
        X4 single 
        X5 single 
        X6 single 
        X7 single 
        U1 single  % (s+c)-dimensional grid of control input 1
        U2 single 
        U3 single 
        X1_next single  % (s+c)-dim grid of next state-1
        X2_next single 
        X3_next single 
        X4_next single 
        X5_next single 
        X6_next single 
        X7_next single 
        
        U1_Opt single  % (s-dim x k) grid of optimal values U* for every stage and states
        U2_Opt single
        U3_Opt single
        %-- ranges for grid generation
        sr_1 % state 1 range vector for grid and interpolant generation 
        sr_2 % "
        sr_3 % "
        sr_4 % "
        sr_5 % "
        sr_6 % "
        sr_7 % "
        cr_1 % control input 1 range vector for grid generation
        cr_2 
        cr_3
        
        h %time step for discrete time system

    end
    
    methods
        %Constructor
        function this = Solver_attitude()
            % states are [w1; w2; w3; q1; q2; q3; q4]
            this.J1 = 2;
            this.J2 = 2.5;
            this.J3 = 3;
            this.dim_U1 = 8;
            this.dim_U2 = 9;
            this.dim_U3 = 10;
            this.U_vector = [-0.01 0 0.01];
        end
        
        %Calculation of optimal matrices
        function obj = run(obj)
            %Preallocation and mesh generation
            size_Xmesh = size(obj.X1_mesh);
            
            obj.U1_Opt = zeros([size_Xmesh,obj.N],'single');
            obj.U1_Opt = zeros([size_Xmesh,obj.N],'single');
            obj.U1_Opt = zeros([size_Xmesh,obj.N],'single');

            % calculate cost to reach next stage
            %obj.J_current_state_fix = g_D(obj);
            calculate_J_current_state_fix(obj);
            % Calculate and store J*NN = h(xi(N)) for all x(N)
            obj.J_next_states_opt = zeros(size(obj.X1),'single');
            %
            calculate_states_next(obj);
            %
            %% obj.stage_J_star = H*X;
            %%
            for k_s=obj.N-1:-1:1
                tic
                calculate_J_U_opt_state_M(obj, k_s);
                fprintf('step %d - %f seconds\n', k_s, toc)
            end
            
        end
        
        function [X1,X2,X3,X4,X5,X6,X7] = a_D_M(obj, X1,X2,X3,X4,X5,X6,X7,U1,U2,U3)
            
            
        end
        
        function calculate_J_current_state_fix(obj)
            obj.J_current_state_fix = obj.Q1*obj.X1.^2 + obj.Q2*obj.X2.^2 + ...
                obj.Q3*obj.X3.^2 + obj.Q4*obj.X4.^2 + ...
                obj.Q5*obj.X5.^2 + obj.Q6*obj.X6.^2 + ...
                obj.Q7*obj.X7.^2 + obj.R1*obj.U1.^2 + ...
                obj.R2*obj.U2.^2 + obj.R3*obj.U3.^2;
        end
        
        function calculate_states_next(obj)
            %calculates next stage (k+1) states
            % X2 = X1 + dt*a_d(X,u) where a_d is the spacecraft dynamics
            
            % Runge-Kutta - 4th order
            % h = dt;
            %             [obj.X1_next, obj.X2_next, obj.X3_next, obj.X4_next,...
            %                 obj.X5_next, obj.X6_next, obj.X7_next] = ...
            %                 spacecraft_dynamics(spacecraft, obj.X1, obj.X2, obj.X3, obj.X4, obj.X5, obj.X6,...
            %                 obj.X7, obj.U1, obj.U2, obj.U3); %X2 = k1
            %
            %             [obj.X1_next, obj.X2_next, obj.X3_next, obj.X4_next,...
            %                 obj.X5_next, obj.X6_next, obj.X7_next] = ...
            %                 spacecraft_dynamics(spacecraft, obj.X1 + obj.X1_next*h/2 , obj.X2 + obj.X2_next*h/2,...
            %                 obj.X3 + obj.X3_next*h/2 , obj.X4 + obj.X4_next*h/2 , obj.X5 + obj.X5_next*h/2 ,...
            %                 obj.X6 + obj.X6_next*h/2 ,...
            %                 obj.X7 + obj.X7_next*h/2 , obj.U1, obj.U2, obj.U3); % k = k2
            %
            %             k = spacecraft_dynamics(spacecraft, (X1 + X2*h/2), U);
            %             X2 = X2 + 2*k; % X2 = k1 + 2*k2
            %             k = spacecraft_dynamics(spacecraft, (X1 + k*h/2), U); %  k = k3
            %             X2 = X2 + 2*k; % X2 = k1 + 2*k2 + 2*k3
            %             k = spacecraft_dynamics(spacecraft, (X1 + k*h), U); %  k = k4
            %             X2 = X2 + k;
            %             X2 = X1 + h/6*(X2);
            
            
            %% taylor expansion
            spacecraft_dynamics_taylor_estimate(obj);
            %
        end
        
        function J = g_D(obj,X1,X2,X3,X4,X5,X6,X7,U1,U2,U3)
            J = obj.Q1*X1.^2 + obj.Q2*X2.^2 + ...
                obj.Q3*X3.^2 + obj.Q4*X4.^2 + ...
                obj.Q5*X5.^2 + obj.Q6*X6.^2 + ...
                obj.Q7*X7.^2 + obj.R1*U1.^2 + ...
                obj.R2*U2.^2 + obj.R3*U3.^2;
        end
        
        function calculate_J_U_opt_state_M(obj, k_s)
            %% CAUTION: this interpolant is only valid for Xmesh
            F = griddedInterpolant(...
                {obj.sr_1, obj.sr_2, obj.sr_3, obj.sr_4,obj.sr_5, obj.sr_6, obj.sr_7}, ...
                obj.J_next_states_opt,'linear');
            
            %find J final for each state and control (X,U) and add it to next state
            %optimum J*
            [obj.J_next_states_opt, U_IDs] = min( obj.J_current_state_fix  + ...
                F(obj.X1,obj.X2,obj.X3,obj.X4,obj.X5,obj.X6,obj.X7) ...
                ,[], obj.dim_U);
            [, , ] = ind2sub(size(obj.J_current_state_fix),id)
        end
        
        function spacecraft_dynamics_taylor_estimate(obj)
            %returns the derivatives x_dot = f(X,u)
            %J is assumed to be diagonal, J12 = J23 = ... = 0
            obj.X1_next = obj.X1 + obj.h*((obj.J2-obj.J3)/obj.J1*obj.X2.*obj.X3 + u1/obj.J1);
            obj.X2_next = obj.X2 + obj.h*((obj.J3-obj.J1)/obj.J2*obj.X3.*obj.X1 + u2/obj.J2);
            obj.X3_next = obj.X3 + obj.h*((obj.J1-obj.J2)/obj.J3*obj.X1.*obj.X2 + u3/obj.J3);
            obj.X4_next = obj.X4 + obj.h*(0.5*(-obj.X1.*obj.X7 -obj.X2.*obj.X6 -obj.X3.*obj.X5));
            obj.X5_next = obj.X5 + obj.h*(0.5*(obj.X2.*obj.X7 -obj.X1.*obj.X6 +obj.X3.*obj.X4));
            obj.X6_next = obj.X6 + obj.h*(0.5*(-obj.X3.*obj.X7 +obj.X1.*obj.X5 +obj.X2.*obj.X4));
            obj.X7_next = obj.X7 + obj.h*(0.5*(obj.X3.*obj.X6 -obj.X2.*obj.X5 +obj.X1.*obj.X4));
        end
        
        function linear_control_response(spacecraft, X0, T_final, dt)
            if nargin < 2
                %sample initial state
                X0 = [1;1;1;...
                    -0.0346073883029131;-0.346079245680143;0.343470774514906;0.872387133925326];
                T_final = 1000;
                dt = 0.01;
            end
            
            
            N = T_final/dt;
            U = [-0.01; 0; 0.01];
            X(:,1) = X0;
            %             qc = [1, 0, 0, 0;...
            %                 0, 1, 0, 0;...
            %                 0, 0, 1, 0;...
            %                 0, 0, 0, 1]; % q command (at origin, is equal to I(4x4) )
            K = [0.5, 0, 0;...
                0, 0.4, 0;
                0, 0, 0.5];
            C = [3, 0, 0;...
                0, 3, 0;
                0, 0, 4];
            %             keyboard;
            tic
            for k_stage=1:N
                % qe = qc*q;
                qe = X(5:7, k_stage);
                q = X(4:7, k_stage);
                w = X(1:3, k_stage);
                U(:,k_stage) = -K*qe - C*w;
                X(:,k_stage+1) = next_stage_states(spacecraft, [w',q'], U(:,k_stage)', dt);
            end
            q_squared_sum = X(4,:).^2 + X(5,:).^2 + X(6,:).^2 + X(7,:).^2; %check quaternions
            %print time and error
            % note: quaternions deviation from (sum(Q.^2) = 1) at T_final is a measure of error in ode solver
            fprintf(...
                'Done - Time elapsed for caculations: %f - states max error: %.1g\n',...
                toc, sqrt((q_squared_sum(end) - 1)))
            
            time_v = linspace(0, T_final, N); %plot time vector
            % plot controls
            figure
            plot(time_v, U(1,:),'--')
            hold on
            plot(time_v, U(2,:),'--')
            plot(time_v, U(3,:),'--')
            grid on
            %plot states
            for n_state = 1:7
                plot(time_v, X(n_state, 1:end-1))
            end
            legend('u1','u2','u3','x1','x2','x3','x4','x5','x6','x7')
            xlabel('time (s)')
        end
        
        function [X1_dot,X2_dot,X3_dot,X4_dot,...
                X5_dot,X6_dot,X7_dot] = spacecraft_dynamics(obj,x1,x2,x3,x4,x5,x6,u1,u2,u3)
            %returns the derivatives x_dot = f(X,u)
            %J is assumed to be diagonal, J12 = J23 = ... = 0
            X1_dot = (obj.J2-obj.J3)/obj.J1*x2.*x3 + u1/obj.J1;
            X2_dot = (obj.J3-obj.J1)/obj.J2*x3.*x1 + u2/obj.J2;
            X3_dot = (obj.J1-obj.J2)/obj.J3*x1.*x2 + u3/obj.J3;
            X4_dot = 0.5*(-x1.*x7 -x2.*x6 -x3.*x5);
            X5_dot = 0.5*(x2.*x7 -x1.*x6 +x3.*x4);
            X6_dot = 0.5*(-x3.*x7 +x1.*x5 +x2.*x4);
            X7_dot = 0.5*(x3.*x6 -x2.*x5 +x1.*x4);
        end
        
        function X_dot = spacecraft_dynamics_list(obj, X,U)
            %returns the derivatives x_dot = f(X,u)
            X_dot = zeros(size(X));
            x1 = X(:,1);
            x2 = X(:,2);
            x3 = X(:,3);
            x4 = X(:,4);
            x5 = X(:,5);
            x6 = X(:,6);
            x7 = X(:,7);
            u1 = U(:,1);
            u2 = U(:,2);
            u3 = U(:,3);
            %J is assumed to be diagonal, J12 = J23 = ... = 0
            X_dot(:,1) = (obj.J2-obj.J3)/obj.J1*x2.*x3 + u1/obj.J1;
            X_dot(:,2) = (obj.J3-obj.J1)/obj.J2*x3.*x1 + u2/obj.J2;
            X_dot(:,3) = (obj.J1-obj.J2)/obj.J3*x1.*x2 + u3/obj.J3;
            X_dot(:,4) = 0.5*(-x1.*x7 -x2.*x6 -x3.*x5);
            X_dot(:,5) = 0.5*(x2.*x7 -x1.*x6 +x3.*x4);
            X_dot(:,6) = 0.5*(-x3.*x7 +x1.*x5 +x2.*x4);
            X_dot(:,7) = 0.5*(x3.*x6 -x2.*x5 +x1.*x4);
        end
        
        function X2 = next_stage_states(spacecraft, X1, U, h)
            %calculates next stage (k+1) states
            % X2 = X1 + dt*a_d(X,u) where a_d is the spacecraft dynamics
            %first order taylor expansion
            %X2 = X1 + dt*spacecraft_dynamics(spacecraft, X1, U);
            
            % Runge-Kutta - 4th order
            % h = dt;
            k1 = spacecraft_dynamics(obj,x1,x2,x3,x4,x5,x6,u1,u2,u3);
            k2 = spacecraft_dynamics(spacecraft, (X1 + k1*h/2), U);
            k3 = spacecraft_dynamics(spacecraft, (X1 + k2*h/2), U);
            k4 = spacecraft_dynamics(spacecraft, (X1 + k3*h), U);
            
            X2 = X1 + h*(k1 + 2*k2 + 2*k3 + k4);
            
        end
        
    end
    
end
