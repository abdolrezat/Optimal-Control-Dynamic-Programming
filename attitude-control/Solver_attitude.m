classdef Solver_attitude < dynamicprops
    %SOLVER_ATTITUDE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        J1
        J2
        J3
        dim_U % dimension number of U, after the 7 state variables
        U_vector % values that are applied as control inputs
        
    end
    
    methods
        %Constructor
        function this = Solver_attitude()
            % states are [w1; w2; w3; q1; q2; q3; q4]
            this.J1 = 2;
            this.J2 = 2.5;
            this.J3 = 3;
            this.dim_U = 8;
            this.U_vector = [-0.01 0 0.01];
        end
        
        %Calculation of optimal matrices
        function obj = run(obj)
            %Preallocation and mesh generation
            size_Xmesh = size(obj.X1_mesh);
            
            obj.stage_J_star = zeros(size_Xmesh,'single');
            obj.u_star = zeros([size_Xmesh,obj.N],'single');
            % Calculate and store J*NN = h(xi(N)) for all x(N)
            %% obj.stage_J_star = H*X;
            %%
            for k_s=obj.N-1:-1:1
                tic
                J_M = J_state_M(obj);
                [obj.stage_J_star, u_star_idx] = min(J_M,[],obj.dim_U);
                % store UMIN in UOPT(N-k,I)
                obj.u_star(:,:,k_s) = obj.U_vector(u_star_idx);
                fprintf('step %d - %f seconds\n', k, toc)
            end
            
        end
        
        function [X1,X2,X3,X4,X5,X6,X7] = a_D_M(obj, X1,X2,X3,X4,X5,X6,X7,U1,U2,U3)
            
            
        end
        function J = g_D(obj,X1,X2,X3,X4,X5,X6,X7,U1,U2,U3)
            J = obj.Q1*X1.^2 + obj.Q2*X2.^2 + ...
                obj.Q3*X3.^2 + obj.Q4*X4.^2 + ...
                obj.Q5*X5.^2 + obj.Q6*X6.^2 + ...
                obj.Q7*X7.^2 + obj.R1*U1.^2 + ...
                obj.R2*U2.^2 + obj.R3*U3.^2;
        end
        
        function J = J_state_M(obj)
            %% CAUSION: this interpolant is only valid for Xmesh
            F = griddedInterpolant(X1,X2,X3,X4,X5,X6,X7, ...
                obj.stage_J_star,'linear');
            % calculate cost to reach next stage
            J = g_D(X1,X2,X3,X4,X5,X6,X7);
            %get next state, replace in Xi to save memory
            [X1, X2, X3, X4, ...
                X5, X6, X7] = a_D_M(X1,X2,X3,X4,X5,X6,X7);
            %find J final for each state and control (X,U) and add it to next state
            %optimum J*
            J = J + F(X1,X2,X3,X4,X5,X6,X7);
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
            
            X2 = X1 + h*(k1 + 2*k2 + 2*k3 + k4)/6;
            
        end
        
        function X2 = next_stage_states_list(spacecraft, X1, U, h)
            %calculates next stage (k+1) states
            % X2 = X1 + dt*a_d(X,u) where a_d is the spacecraft dynamics
            %first order taylor expansion
            %X2 = X1 + dt*spacecraft_dynamics(spacecraft, X1, U);
            
            % Runge-Kutta - 4th order
            % h = dt;
            k1 = spacecraft_dynamics_list(spacecraft, X1, U);
            k2 = spacecraft_dynamics_list(spacecraft, (X1 + k1*h/2), U);
            k3 = spacecraft_dynamics_list(spacecraft, (X1 + k2*h/2), U);
            k4 = spacecraft_dynamics_list(spacecraft, (X1 + k3*h), U);
            
            X2 = X1 + h*(k1 + 2*k2 + 2*k3 + k4)/6;
            
        end
    end
    
end

