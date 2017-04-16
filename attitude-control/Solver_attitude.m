classdef Solver_attitude < dynamicprops
    %SOLVER_ATTITUDE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        J1
        J2
        J3
    end
    
    methods
        %Constructor
        function this = Solver_attitude()
            % states are [w1; w2; w3; q1; q2; q3; q4]
            this.J1 = 2;
            this.J2 = 2.5;
            this.J3 = 3;
        end
        
        %Calculation of optimal matrices
        function obj = run(obj)

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
        
        function X_dot = spacecraft_dynamics(obj, X,U)
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
            k1 = spacecraft_dynamics(spacecraft, X1, U);
            k2 = spacecraft_dynamics(spacecraft, (X1 + k1*h/2), U); 
            k3 = spacecraft_dynamics(spacecraft, (X1 + k2*h/2), U);
            k4 = spacecraft_dynamics(spacecraft, (X1 + k3*h), U);
            
            X2 = X1 + h*(k1 + 2*k2 + 2*k3 + k4)/6;
            
        end
    end
    
end

