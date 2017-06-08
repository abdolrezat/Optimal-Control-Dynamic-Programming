classdef Solver_pos_att < handle
    %SOLVER_POSITION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        N % number of stages
        Mass % Mass
        InertiaM % Moment of Inertia Matrix
        J1 % element 1 of Inertia Matrix
        J2 % ...
        J3 % ..
        
        size_state_mat
        
        T_vector % Thruster Forces that are applied as control inputs
        
        %pos
        v_min
        v_max
        n_mesh_v
        
        x_min
        x_max
        n_mesh_x
        
        %att
        w_min
        w_max
        n_mesh_w
        
        theta1_min
        theta1_max
        theta2_min
        theta2_max
        theta3_min
        theta3_max
        n_mesh_t
        
        %
        Qx1
        Qx2
        Qx3
        Qv1
        Qv2
        Qv3
        
        Qt1
        Qt2
        Qt3
        Qw1
        Qw2
        Qw3
        
        R1
        R2
        R3
        
        T_final
        h
        N_stage
        defaultX0
        
        %Thruster and Torque forces
        F_Thr0 %Thruster 0 - Force
        F_Thr1
        F_Thr2
        F_Thr3
        F_Thr4
        F_Thr5
        F_Thr6
        F_Thr7
        F_Thr8
        F_Thr9
        F_Thr10
        F_Thr11
        
        T_dist %Thruster placement distanceS
        %
        %Optimum values for each Thruster
        Opt_F_Thr0
        Opt_F_Thr1
        Opt_F_Thr2
        Opt_F_Thr3
        Opt_F_Thr4
        Opt_F_Thr5
        Opt_F_Thr6
        Opt_F_Thr7
        Opt_F_Thr8
        Opt_F_Thr9
        Opt_F_Thr10
        Opt_F_Thr11
        
    end
    
    methods
        function this = Solver_pos_att()
            if nargin < 1
                
                %pos
                this.v_min = -0.5;
                this.v_max = +0.5;
                this.n_mesh_v = 50;
                
                this.x_min = -0.5;
                this.x_max = 0.5;
                this.n_mesh_x = 50;
                
                %att
                this.w_min = -deg2rad(50);
                this.w_max = -deg2rad(-50);
                this.n_mesh_w = 5;
                
                this.theta1_min = -5; %angles of rotation about y-axis (pitch)
                this.theta1_max = 5;
                this.theta2_min = -6; %angles of rotation about z-axis (yaw)
                this.theta2_max = 6;
                this.theta3_min = -7; %x-axis rotation
                this.theta3_max = 7;
                this.n_mesh_t = 5;
                %
                
                this.Mass = 4.16;
                inertia(1,1) =  0.02836 + 0.00016;
                inertia(2,1) =  0.026817 + 0.00150;
                inertia(3,1) =  0.023 + 0.00150;
                inertia(4,1) = -0.0000837;
                inertia(5,1) =  0.000014;
                inertia(6,1) = -0.00029;
                this.InertiaM = [inertia(1,1)  inertia(4,1)  inertia(5,1);...
                    inertia(4,1)  inertia(2,1)  inertia(6,1);...
                    inertia(5,1)  inertia(6,1)  inertia(3,1)];
                
                this.J1 = this.InertiaM(1);
                this.J2 = this.InertiaM(5);
                this.J3 = this.InertiaM(9);
                
                this.Qx1 = 6;
                this.Qx2 = 6;
                this.Qx3 = 6;
                this.Qv1 = 6;
                this.Qv2 = 6;
                this.Qv3 = 6;
                
                this.Qt1 = .5;
                this.Qt2 = .5;
                this.Qt3 = .5;
                this.Qw1 = .5;
                this.Qw2 = .5;
                this.Qw3 = .5;
                
                this.R1 = 0.1;
                this.R2 = 0.1;
                this.R3 = 0.1;
                
                this.T_final = 10;
                this.h = 0.005;
            end
            
            this.N_stage = this.T_final/this.h;
            
            if(~isinteger( this.N_stage))
                this.N_stage = ceil(this.N_stage);
                this.T_final = this.h*this.N_stage;
                warning('T_final is not a factor of h (dt), increasing T_final to %.2f\n',this.T_final)
            end
            
            q0 = [];
            this.defaultX0 = [0;0;0;0;0;0;q0;0;0;0]; % [x;v;q;w]
            
            %Thruster Forces
            Thruster_max_F = 0.13; % (N)
            this.T_dist = 9.65E-2; % (meters)
            %x-direction
            this.F_Thr0 = [0 Thruster_max_F];
            this.F_Thr1 = [0 Thruster_max_F];
            
            this.F_Thr6 = -[0 Thruster_max_F];
            this.F_Thr7 = -[0 Thruster_max_F];
            
            %y-direction
            this.F_Thr2 = [0 Thruster_max_F];
            this.F_Thr3 = [0 Thruster_max_F];
            
            this.F_Thr8 = -[0 Thruster_max_F];
            this.F_Thr9 = -[0 Thruster_max_F];
            
            %z-direction
            this.F_Thr4 = [0 Thruster_max_F];
            this.F_Thr5 = [0 Thruster_max_F];
            
            this.F_Thr10 = -[0 Thruster_max_F];
            this.F_Thr11 = -[0 Thruster_max_F];
            
            
        end
        
        function simplified_run(obj)
            
            %% mesh generation
            s_x1 = sym_linspace(obj, obj.x_min, obj.x_max, obj.n_mesh_x);
            s_x2 = sym_linspace(obj, obj.x_min, obj.x_max, obj.n_mesh_x);
            s_x3 = sym_linspace(obj, obj.x_min, obj.x_max, obj.n_mesh_x);
            
            s_v1 = sym_linspace(obj, obj.v_min, obj.v_max, obj.n_mesh_v);
            s_v2 = sym_linspace(obj, obj.v_min, obj.v_max, obj.n_mesh_v);
            s_v3 = sym_linspace(obj, obj.v_min, obj.v_max, obj.n_mesh_v);
            
            s_t1 = sym_linspace(obj, deg2rad(obj.theta1_min), deg2rad(obj.theta1_max), obj.n_mesh_t);
            s_t2 = sym_linspace(obj, deg2rad(obj.theta2_min), deg2rad(obj.theta2_max), obj.n_mesh_t);
            s_t3 = sym_linspace(obj, deg2rad(obj.theta3_min), deg2rad(obj.theta3_max), obj.n_mesh_t);
            
            s_w1 = sym_linspace(obj, obj.w_min, obj.w_max, obj.n_mesh_w);
            s_w2 = sym_linspace(obj, obj.w_min, obj.w_max, obj.n_mesh_w);
            s_w3 = sym_linspace(obj, obj.w_min, obj.w_max, obj.n_mesh_w);
            
            
            %% initialization
            %calculating J fixed
            [f0_allcomb,f1_allcomb,f6_allcomb,f7_allcomb] = ...
                get_allcomb_vectors(obj, obj.F_Thr0, obj.F_Thr1, obj.F_Thr6, obj.F_Thr7);

            J_current_1 = J_current_reshaped(obj, s_x1,s_v1,s_t1,s_w1,...
                f0_allcomb,  f1_allcomb,  f6_allcomb,  f7_allcomb,...
                obj.Qx1,obj.Qv1,obj.Qt1,obj.Qw1,obj.R1);
            F1 = griddedInterpolant({s_x1,s_v1,s_t1,s_w1},...
                zeros(obj.n_mesh_x,obj.n_mesh_v,obj.n_mesh_t,obj.n_mesh_w,'single'),'linear');
            
%             J_current_2 = J_current_reshaped(obj, s_x2,s_v2,s_t2,s_w2,...
%                 obj.F_Thr2, obj.F_Thr3, obj.F_Thr8, obj.F_Thr9,...
%                 obj.Qx2,obj.Qv2,obj.Qt2,obj.Qw2,obj.R2);
%             F2 = griddedInterpolant({s_x2,s_v2,s_t2,s_w2},...
%                 zeros(obj.n_mesh_x,obj.n_mesh_v,obj.n_mesh_t,obj.n_mesh_w,'single'),'linear');
%             
%             J_current_3 = J_current_reshaped(obj, s_x3,s_v3,s_t3,s_w3,...
%                 obj.F_Thr4, obj.F_Thr5, obj.F_Thr10, obj.F_Thr11,...
%                 obj.Qx3,obj.Qv3,obj.Qt3,obj.Qw3,obj.R3);
%             F3 = griddedInterpolant({s_x3,s_v3,s_t3,s_w3},...
%                 zeros(obj.n_mesh_x,obj.n_mesh_v,obj.n_mesh_t,obj.n_mesh_w,'single'),'linear');
            
            %calculating next stage states
            %
            [x1_next,v1_next,t1_next,w1_next] = next_stage_states_simplified(obj, s_x1,s_v1,s_t1,s_w1,...
                f0_allcomb, f1_allcomb, f6_allcomb, f7_allcomb, obj.J1);
%             [x2_next,v2_next,t2_next,w2_next] = next_stage_states_simplified(obj, s_x2,s_v2,s_t2,s_w2,...
%                 obj.F_Thr2, obj.F_Thr3, obj.F_Thr8, obj.F_Thr9, obj.J2);
%             [x3_next,v3_next,t3_next,w3_next] = next_stage_states_simplified(obj, s_x3,s_v3,s_t3,s_w3,...
%                 obj.F_Thr4, obj.F_Thr5, obj.F_Thr10, obj.F_Thr11, obj.J3);
            
            %beginning (reverse) stage calculations
%             U_Optimal_id = zeros(obj.n_mesh_x,obj.n_mesh_v,obj.n_mesh_t,obj.n_mesh_w,'single');
            
            whandle = waitbar(0,'Calculation in Progress...');
            for k_s = obj.N_stage-1:-1:1
                tic
                [F1.Values, U_Optimal_id] = min( J_current_1 + F1(x1_next,v1_next,t1_next,w1_next), [], 5);

%                 
%                 [val, obj.Opt_F_Thr9] = min( J_current_2 + F2(x2_next,v2_next,t2_next,w2_next), [], 8);
%                 [val, obj.Opt_F_Thr8] = min(  val, [], 7);
%                 [val, obj.Opt_F_Thr3] = min(  val, [], 6);
%                 [F2.Values, obj.Opt_F_Thr2] = min( val, [], 5);
%                 
%                 [val, obj.Opt_F_Thr11] = min( J_current_3 + F3(x3_next,v3_next,t3_next,w3_next), [], 8);
%                 [val, obj.Opt_F_Thr10] = min(  val, [], 7);
%                 [val, obj.Opt_F_Thr5] = min(  val, [], 6);
%                 [F3.Values, obj.Opt_F_Thr4] = min( val, [], 5);
                
                waitbar( 1 - k_s/obj.N_stage, whandle);
                fprintf('step %d - %f seconds\n', k_s, toc)
            end
            
%             disp(unique(U_Optimal_id)')
            %get U* Optimal idx
            obj.Opt_F_Thr0 = f0_allcomb(U_Optimal_id);
            obj.Opt_F_Thr1 = f1_allcomb(U_Optimal_id);
            obj.Opt_F_Thr6 = f6_allcomb(U_Optimal_id);
            obj.Opt_F_Thr7 = f7_allcomb(U_Optimal_id);
            
            %set U* Optimal values
            obj.Opt_F_Thr0 = griddedInterpolant({s_x1,s_v1,s_t1,s_w1},...
                obj.Opt_F_Thr0,...
                'nearest');
            obj.Opt_F_Thr1 = griddedInterpolant({s_x1,s_v1,s_t1,s_w1},...
                obj.Opt_F_Thr1,...
                'nearest');
            obj.Opt_F_Thr6 = griddedInterpolant({s_x1,s_v1,s_t1,s_w1},...
                obj.Opt_F_Thr6,...
                'nearest');
            obj.Opt_F_Thr7 = griddedInterpolant({s_x1,s_v1,s_t1,s_w1},...
                obj.Opt_F_Thr7,...
                'nearest');
            
            %copy one channel to two others (TEST)
            
            
            obj.Opt_F_Thr9 = obj.Opt_F_Thr7;
            obj.Opt_F_Thr8 = obj.Opt_F_Thr6;
            obj.Opt_F_Thr3 = obj.Opt_F_Thr1;
            obj.Opt_F_Thr2 = obj.Opt_F_Thr0;
            
            obj.Opt_F_Thr11 = obj.Opt_F_Thr7;
            obj.Opt_F_Thr10 = obj.Opt_F_Thr6;
            obj.Opt_F_Thr5 = obj.Opt_F_Thr1;
            obj.Opt_F_Thr4 = obj.Opt_F_Thr0;
%             

%             
%             obj.Opt_F_Thr2 = griddedInterpolant({s_x2,s_v2,s_t2,s_w2},...
%                 obj.F_Thr2(obj.Opt_F_Thr2),...
%                 'nearest');
%             obj.Opt_F_Thr3 = griddedInterpolant({s_x2,s_v2,s_t2,s_w2},...
%                 obj.F_Thr3(obj.Opt_F_Thr3),...
%                 'nearest');
%             obj.Opt_F_Thr8 = griddedInterpolant({s_x2,s_v2,s_t2,s_w2},...
%                 obj.F_Thr8(obj.Opt_F_Thr8),...
%                 'nearest');
%             obj.Opt_F_Thr9 = griddedInterpolant({s_x2,s_v2,s_t2,s_w2},...
%                 obj.F_Thr9(obj.Opt_F_Thr9),...
%                 'nearest');
%             
%             obj.Opt_F_Thr4 = griddedInterpolant({s_x3,s_v3,s_t3,s_w3},...
%                 obj.F_Thr4(obj.Opt_F_Thr4),...
%                 'nearest');
%             obj.Opt_F_Thr5 = griddedInterpolant({s_x3,s_v3,s_t3,s_w3},...
%                 obj.F_Thr5(obj.Opt_F_Thr5),...
%                 'nearest');
%             obj.Opt_F_Thr10 = griddedInterpolant({s_x3,s_v3,s_t3,s_w3},...
%                 obj.F_Thr10(obj.Opt_F_Thr10),...
%                 'nearest');
%             obj.Opt_F_Thr11 = griddedInterpolant({s_x3,s_v3,s_t3,s_w3},...
%                 obj.F_Thr11(obj.Opt_F_Thr11),...
%                 'nearest');
%             
            %finish up
            close(whandle)
            fprintf('stage calculations complete.\n')
            
        end
        
        function [x_next,v_next,t_next,w_next] = next_stage_states_simplified(obj, X, V, T, W, f1,f2,f6,f7 ,J)
            % reshape
            X = reshape(X,[obj.n_mesh_x 1]  );
            V = reshape(V,[1 obj.n_mesh_v]  );
            T = reshape(T,[1 1 obj.n_mesh_t]  );
            W = reshape(W,[1 1 1 obj.n_mesh_w]  );
            f1 = reshape(f1,[1 1 1 1 length(f1)]  );
            f2 = reshape(f2,[1 1 1 1 length(f2)]  );
            f6 = reshape(f6,[1 1 1 1 length(f6)]  );
            f7 = reshape(f7,[1 1 1 1 length(f7)]  );
            
            % ODE solve
            x_next = RK4_x(obj, X, V);
            v_next = RK4_v(obj, V, f1,f2,f6,f7);
            t_next = RK4_t(obj, T, W);
            w_next = RK4_w(obj, W, f1,f2,f6,f7 , J);
            

            %repmat each matrix to full size, as required for F inputs
            x_next = repmat(x_next,[1 1 obj.n_mesh_t obj.n_mesh_w length(f1)]);
            v_next = repmat(v_next,[obj.n_mesh_x 1 obj.n_mesh_t obj.n_mesh_w 1]);
            t_next = repmat(t_next,[obj.n_mesh_x obj.n_mesh_v 1 1 length(f1)]);
            w_next = repmat(w_next,[obj.n_mesh_x obj.n_mesh_v obj.n_mesh_t 1 1]);    
            end
        
        function X2 = RK4_x(obj, X1, V)
            % Runge-Kutta - 4th order
            % h = dt;
            k1 = xdynamics(obj, V);
%             k2 = xdynamics(obj,(V + k1*obj.h/2));
%             k3 = xdynamics(obj,(V + k2*obj.h/2));
%             k4 = xdynamics(obj,(V + k3*obj.h));
            
            X2 = X1 + obj.h*k1;
            
        end
        
        function x_dot = xdynamics(~,v)
            x_dot = v;
        end
        
        function V2 = RK4_v(obj, V1, f1,f2,f6,f7) % does not need RK4, ki's are equal
            % Runge-Kutta - 4th order
            % h = dt;
            k1 = vdynamics(obj, V1 , f1,f2,f6,f7);
%             k2 = vdynamics(obj,(V1 + k1*obj.h/2), f1,f2,f6,f7);
%             k3 = vdynamics(obj,(V1 + k2*obj.h/2), f1,f2,f6,f7);
%             k4 = vdynamics(obj,(V1 + k3*obj.h), f1,f2,f6,f7);
            
            V2 = V1 + obj.h*k1;
        end
        
        function v_dot = vdynamics(obj, ~, f1,f2,f6,f7)
            v_dot = (f1+f2+f6+f7)/obj.Mass;
        end
        
        function T2 = RK4_t(obj, T1, W1)
            %calculates next stage (k+1) states
            % X2 = X1 + dt*a_d(X,u) where a_d is the spacecraft dynamics
            %first order taylor expansion
            %X2 = X1 + dt*spacecraft_dynamics(spacecraft, X1, U);
            
            % Runge-Kutta - 4th order
            % h = dt;
            k1 = tdynamics(obj,W1);
%             k2 = tdynamics(obj,(W1 + k1*obj.h/2));
%             k3 = tdynamics(obj,(W1 + k2*obj.h/2));
%             k4 = tdynamics(obj,(W1 + k3*obj.h));
            
            T2 = T1 + obj.h*k1;
        end
        
        
        function t_dot = tdynamics(~,w)
            t_dot = w;
        end
        
        function W2 = RK4_w(obj, w, f1,f2,f6,f7, J)
            %calculates next stage (k+1) states
            % X2 = X1 + dt*a_d(X,u) where a_d is the spacecraft dynamics
            %first order taylor expansion
            %X2 = X1 + dt*spacecraft_dynamics(spacecraft, X1, U);
            
            % Runge-Kutta - 4th order
            % h = dt;
            k1 = wdynamics(obj,w , f1,f2,f6,f7, J);
%             k2 = wdynamics(obj,(w + k1*obj.h/2), f1,f2,f6,f7, J);
%             k3 = wdynamics(obj,(w + k2*obj.h/2), f1,f2,f6,f7, J);
%             k4 = wdynamics(obj,(w + k3*obj.h), f1,f2,f6,f7, J);
            
            W2 = w + obj.h*(k1);
        end
        
        
        function w_dot = wdynamics(obj,~, f1,f2,f6,f7, J)
            w_dot = (f1*obj.T_dist +f2*(-obj.T_dist) + ...
                f6*obj.T_dist +f7*(-obj.T_dist) )/J;
        end
        
        function [f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11] = get_thruster_on_off_optimal(obj,x,v,t,w,R0,V0,q)
            % gets the optimal on/off state of thrusters, method run() or
            % simplified_run() must be run before calling this function,
            % the inputs x,v are position and velocity of object B
            % to target A in RSW frame
            
            %transform vectors to body frame of reference
            rotM_RSW2ECI = RSW2ECI(obj, R0, V0);
            rotM_ECI2body = ECI2body(obj,q);
            
            x = rotM_ECI2body*(rotM_RSW2ECI* x');
            v = rotM_ECI2body*(rotM_RSW2ECI* v');
            
            %
            t_x = t(1); %rotation about x_axis
            t_y = t(2);
            t_z = t(3);
            
            w_x = w(1); %rotational speed about x_axis
            w_y = w(2);
            w_z = w(3);
            
            x1 = x(1);
            x2 = x(2);
            x3 = x(3);
            
            v1 = v(1);
            v2 = v(2);
            v3 = v(3);
            
            f0 = obj.Opt_F_Thr0(x1,v1, t_y, w_y );
            f1 = obj.Opt_F_Thr1(x1,v1, t_y, w_y );
            f6 = obj.Opt_F_Thr6(x1,v1, t_y, w_y );
            f7 = obj.Opt_F_Thr7(x1,v1, t_y, w_y );
            
            f2 = obj.Opt_F_Thr2(x2,v2, t_z, w_z );
            f3 = obj.Opt_F_Thr3(x2,v2, t_z, w_z );
            f8 = obj.Opt_F_Thr8(x2,v2, t_z, w_z );
            f9 = obj.Opt_F_Thr9(x2,v2, t_z, w_z );
            
            f4 = obj.Opt_F_Thr4(x3,v3, t_x, w_x );
            f5 = obj.Opt_F_Thr5(x3,v3, t_x, w_x );
            f10 = obj.Opt_F_Thr10(x3,v3, t_x, w_x );
            f11 = obj.Opt_F_Thr11(x3,v3, t_x, w_x );
            
        end
        
        
        function get_optimal_path(obj)
            global mu
            mu = 398600;
            if nargin < 2
%                 X0 = obj.defaultX0;
                %   Prescribed initial state vector of chaser B in the co-moving frame:
                dr0 = [-1  0  0];
                dv0 = [0 0 0];
                q0 = [0 0 0 1];
                w0 = [0 0 0];
                X0 = [dr0 dv0 q0 w0]';
                tf  = obj.T_final;
                N_total_sim = obj.N_stage;
            end
            
            tspan = 0:obj.h:tf;
            X_ode45 = zeros(N_total_sim, 13);
            F_Th_Opt = zeros(N_total_sim, 12);
            Force_Moment_log = zeros(N_total_sim, 6);
            X_ode45(1,:) = X0;
            % Calculate the target initial state vector
            [R0,V0] = get_target_R0V0(obj);
            tic
            for k_stage=1:N_total_sim-1
                %determine F_Opt each Thruster
                X_stage = X_ode45(k_stage,:);
                
                x_stage = X_stage(1:3);
                v_stage = X_stage(4:6);
                t_stage = [2*asin(X_stage(9));... %angle about x-axis
                    2*asin(X_stage(8));... %y-axis
                    2*asin(X_stage(7))]; %z-axis
                w_stage = X_stage(11:13);
                q_stage = X_stage(7:10);
                [f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11] = ...
                    get_thruster_on_off_optimal(obj, x_stage, v_stage, t_stage, w_stage, R0,V0,q_stage);

                % calculate moments (U_M) and directional forces (a_* |x,y,z|)
                [U_M, a_x, a_y, a_z] = to_Moments_Forces(obj,f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,R0,V0,q_stage);
                %log
                F_Th_Opt(k_stage,:) = [f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11];
                Force_Moment_log(k_stage,:) = [a_x, a_y, a_z, U_M'];
                %
                [~,X_temp] = ode45(@ode_eq,[tspan(k_stage), tspan(k_stage+1)], X_stage);
                X_ode45(k_stage+1,:) = X_temp(end,:);
            end
        toc
        T_ode45 = tspan(1:end-1)';
        %plot Thruster Firings
        ylim_thr = [-.15 .15];
        figure;
        title('Thruster Firings')
        subplot(4,3,1)
        plot(T_ode45, F_Th_Opt(:,1))
        title('#0 (x)')
        grid on
        ylim(ylim_thr)

        subplot(4,3,2)
        plot(T_ode45, F_Th_Opt(:,3))
        title('#2 (y)')    
        grid on
        ylim(ylim_thr)
        
        subplot(4,3,3)
        plot(T_ode45, F_Th_Opt(:,5))
        title('#4 (z)')
        grid on
        ylim(ylim_thr)
        
        subplot(4,3,4)
        plot(T_ode45, F_Th_Opt(:,2))
        title('#1 (x)')
        grid on
        ylim(ylim_thr)
        
        subplot(4,3,5)
        plot(T_ode45, F_Th_Opt(:,4))
        title('#3 (y)')
        grid on
        ylim(ylim_thr)
        
        subplot(4,3,6)
        plot(T_ode45, F_Th_Opt(:,6))
        title('#5 (z)')
        grid on
        ylim(ylim_thr)
        
        subplot(4,3,7)
        plot(T_ode45, F_Th_Opt(:,7))
        title('#6 (-x)')
        grid on
        ylim(ylim_thr)
        
        subplot(4,3,8)
        plot(T_ode45, F_Th_Opt(:,9))
        title('#8 (-y)')
        grid on
        ylim(ylim_thr)
        
        subplot(4,3,9)
        plot(T_ode45, F_Th_Opt(:,11))
        title('#10 (-z)')
        grid on
        ylim(ylim_thr)
        
        subplot(4,3,10)
        plot(T_ode45, F_Th_Opt(:,8))
        title('#7 (-x)')
        grid on
        ylim(ylim_thr)
        
        subplot(4,3,11)
        plot(T_ode45, F_Th_Opt(:,10))
        title('#9 (-y)')
        grid on
        ylim(ylim_thr)
        
        subplot(4,3,12)
        plot(T_ode45, F_Th_Opt(:,12))
        title('#11 (-z)')
        grid on
        ylim(ylim_thr)
        
        %plot Moments
        ylim_F = [-0.3 0.3]/obj.Mass;
        ylim_M = [-0.3 0.3]*obj.T_dist;
        figure;
        title('Forces and Moments')
        
        subplot(3,2,1)
        plot(T_ode45, Force_Moment_log(:,1))
        title('acceleration x-direction')
        grid on
        ylim(ylim_F)
        
        subplot(3,2,3)
        plot(T_ode45, Force_Moment_log(:,2))
        title('acceleration y-direction')
        grid on
        ylim(ylim_F)
        
        subplot(3,2,5)
        plot(T_ode45, Force_Moment_log(:,3))
        title('acceleration z-direction')
        grid on
        ylim(ylim_F)
        
        subplot(3,2,2)
        plot(T_ode45, Force_Moment_log(:,4))
        title('Moment x-direction')
        grid on
        ylim(ylim_M)
        
        subplot(3,2,4)
        plot(T_ode45, Force_Moment_log(:,5))
        title('Moment y-direction')
        grid on
        ylim(ylim_M)
        
        subplot(3,2,6)
        plot(T_ode45, Force_Moment_log(:,6))
        title('Moment z-direction')
        grid on
        ylim(ylim_M)
        
        % plot states - pos
        figure;
        title('states - position')
        plot(T_ode45, X_ode45(:,1))
        hold on
        plot(T_ode45, X_ode45(:,2))
        plot(T_ode45, X_ode45(:,3))
        grid on
        legend('x1','x2','x3')
        
        % plot states - v
        figure;
        title('states - velocity')
        plot(T_ode45, X_ode45(:,4))
        hold on
        plot(T_ode45, X_ode45(:,5))
        plot(T_ode45, X_ode45(:,6))
        grid on
        legend('v1','v2','v3')
        
        % plot states - q
        figure;
        title('states - quaternions')
        plot(T_ode45, X_ode45(:,7))
        hold on
        plot(T_ode45, X_ode45(:,8))
        plot(T_ode45, X_ode45(:,9))
        plot(T_ode45, X_ode45(:,10))
        grid on
        legend('q1','q2','q3','q4')
        
        % plot states - q
        figure;
        title('states - rotational speeds')
        plot(T_ode45, X_ode45(:,11))
        hold on
        plot(T_ode45, X_ode45(:,12))
        plot(T_ode45, X_ode45(:,13))
        grid on
        legend('w1','w2','w3')
        
        %function declarations
            function x_dot = ode_eq(t,X1)
                x_dot = system_dynamics(t,X1);
                x_dot = x_dot';
                function X_dot = system_dynamics(t,X)
                    x1 = X(1);
                    x2 = X(2);
                    x3 = X(3);
                    v1 = X(4);
                    v2 = X(5);
                    v3 = X(6);
                    q1 = X(7);
                    q2 = X(8);
                    q3 = X(9);
                    q4 = X(10);
                    w1 = X(11);
                    w2 = X(12);
                    w3 = X(13);
                    w_vector = X(11:13);
                    %--- differential equations -------------------------
                    % pre-computations
                    [R,V] = update_RV_target(obj, R0, V0, t);
                    norm_R = (R*R')^.5; %norm R
                    RdotV = sum(R.*V); %dot product
                    crossRV = [R(2).*V(3)-R(3).*V(2); % cross product of R and V
                        R(3).*V(1)-R(1).*V(3);
                        R(1).*V(2)-R(2).*V(1)];
                    H  = (crossRV'*crossRV)^.5 ; %norm(crossRV);
                    
                    % CW-equations
                    % position - x
                    X_dot(1) = v1;
                    X_dot(2) = v2;
                    X_dot(3) = v3;
                    
                    % position - v (a_x,y,z are in RSW frame of reference)
                    X_dot(4) =  (2*mu/norm_R^3 + H^2/norm_R^4)*x1 - 2*RdotV/norm_R^4*H*x2 + 2*H/norm_R^2*v2 ...
                        + a_x;
                    X_dot(5) = -(mu/norm_R^3 - H^2/norm_R^4)*x2 + 2*RdotV/norm_R^4*H*x1 - 2*H/norm_R^2*v1 ...
                        + a_y;
                    X_dot(6) = -mu/norm_R^3*x3 ...
                        + a_z;
                    
                    % attitude - q
                    X_dot(7) = 0.5*(w3.*q2 -w2.*q3 +w1.*q4);
                    X_dot(8) = 0.5*(-w3.*q1 +w1.*q3 +w2.*q4);
                    X_dot(9) = 0.5*(w2.*q1 -w1.*q2 +w3.*q4);
                    X_dot(10) = 0.5*(-w1.*q1 -w2.*q2 -w3.*q3);
                    
                    % attitude - w
                    w_dot = obj.InertiaM\(U_M - cross(w_vector, obj.InertiaM*w_vector));
                    X_dot(11) = w_dot(1);
                    X_dot(12) = w_dot(2);
                    X_dot(13) = w_dot(3);
                end
            end
        end
        
   
        
        function [R0,V0] = get_target_R0V0(obj)
            global mu
            RE  = 6378;
            %...Input data:
            %   Prescribed initial orbital parameters of target A:
            rp    = RE + 300;
            e     = 0.1;
            i     = 0;
            RA    = 0;
            omega = 0;
            theta = 0;
            %   Additional computed parameters:
            ra = rp*(1 + e)/(1 - e);
            h_  = sqrt(2*mu*rp*ra/(ra + rp));
            a  = (rp + ra)/2;
            T  = 2*pi/sqrt(mu)*a^1.5;
            n  = 2*pi/T;
            [R0,V0] = sv_from_coe([h_ e RA i omega theta],mu);
        end
        
        function [R2,V2] = update_RV_target(~,R0,V0,t)
            %Updates the state vectors of the target sat
            % mu - gravitational parameter (km^3/s^2)
            % R0 - initial position vector (km)
            % V0 - initial velocity vector (km/s)
            % t - elapsed time (s)
            % R - final position vector (km)
            % V - final velocity vector (km/s)
            global mu
            %...Magnitudes of R0 and V0:
            r0 = norm(R0);
            v0 = norm(V0);
            %...Initial radial velocity:
            vr0 = dot(R0, V0)/r0;
            %...Reciprocal of the semimajor axis (from the energy equation):
            alpha = 2/r0 - v0^2/mu;
            %...Compute the universal anomaly:
            x = kepler_U(t, r0, vr0, alpha);
            %...Compute the f and g functions:
            [f, g] = f_and_g(x, t, r0, alpha);
            %...Compute the final position vector:
            R2 = f*R0 + g*V0;
            %...Compute the magnitude of R:
            r2 = norm(R2);
            %...Compute the derivatives of f and g:
            [fdot, gdot] = fDot_and_gDot(x, r2, r0, alpha);
            %...Compute the final velocity:
            V2 = fdot*R0 + gdot*V0;
        end
        
        function J_current_M = J_current_reshaped(obj,x,v,t,w,f1,f2,f3,f4,Qx,Qv,Qt,Qw,R)
            x = reshape(x,[obj.n_mesh_x 1]  );
            v = reshape(v,[1 obj.n_mesh_v]  );
            t = reshape(t,[1 1 obj.n_mesh_t]  );
            w = reshape(w,[1 1 1 obj.n_mesh_w]  );
            f1 = reshape(f1,[1 1 1 1 length(f1)]  );
            f2 = reshape(f2,[1 1 1 1 length(f2)]  );
            f3 = reshape(f3,[1 1 1 1 length(f3)]  );
            f4 = reshape(f4,[1 1 1 1 length(f4)]  );
%             
            J_current_M = single(Qx * x.^2 + Qv * v.^2 + Qw * w.^2 + Qt * t.^2 +...
                ( R* f1.^2 +  R* f2.^2 +  R* f3.^2 +  R* f4.^2)  );
        end
        
        function [U_M, a_x, a_y, a_z] = to_Moments_Forces(obj,f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,R0,V0,q)
            % Moments
            U_M_y = (f0-f1+f6-f7)*obj.T_dist;
            U_M_z = (f2-f3+f8-f9)*obj.T_dist;
            U_M_x = (f4-f5+f10-f11)*obj.T_dist;
            U_M = [U_M_x; U_M_y; U_M_z];
            
            % Forces (expressed in body frame of reference)
            a_x_body = (f0+f1+f6+f7)/obj.Mass;
            a_y_body = (f2+f3+f8+f9)/obj.Mass;
            a_z_body = (f4+f5+f10+f11)/obj.Mass;
            % transform vectors
            rotM_RSW2ECI = RSW2ECI(obj, R0, V0);
            rotM_ECI2body = ECI2body(obj,q);
            
            accM =  rotM_RSW2ECI\(rotM_ECI2body\[a_x_body a_y_body a_z_body]');
            a_x = accM(1);
            a_y = accM(2);
            a_z = accM(3);

        end
        
        function qrotMat = ECI2body(~, q)
            qrotMat = [1-2*(q(2)^2+q(3)^2), 2*(q(1)*q(2)+q(3)*q(4)), 2*(q(1)*q(3)-q(2)*q(4));...
                2*(q(2)*q(1)-q(3)*q(4)), 1-2*(q(1)^2+q(3)^2), 2*(q(2)*q(3) +q(1)*q(4));...
                2*(q(3)*q(1)+q(2)*q(4)), 2*(q(3)*q(2)-q(1)*q(4)), 1-2*(q(1)^2+q(2)^2)];
        end
        
        function rotMat = RSW2ECI(~, pos, vel)
            % rotMat = RSW2ECI(pos, vel);
            % Creates a rotation matrix which transforms RSW vectors to ECI vectors.
            % ECIvec = rotMat  *   RSW;
            % 3x1    = 3x3         3x1;
            % Inputs:
            %   pos:   ECI position vector
            %   vel:   ECI velocity vector
            % Outpus:
            %   rotMat: 3x3 rotation matrix from RSW to ECI
            
            R = pos/norm(pos);
            W = cross(pos,vel)/norm(cross(pos,vel));
            S = cross(W,R);
            
            rotMat = [R' S' W'];
        end
        
        function  [f1_a,f2_a,f3_a,f4_a] = get_allcomb_vectors(~,f1,f2,f3,f4)
            [f1,f2,f3,f4] = ndgrid(f1,f2,f3,f4);
             f1_a = f1(:); 
             f2_a = f2(:); 
             f3_a = f3(:); 
             f4_a = f4(:); 
        end
        
        function v = sym_linspace(~,a,b,n)
            if(a>0)
                error('minimum states are not negative, use normal linspace')
            end
            if(mod(n,2) == 0)
                v_1 = linspace(a,0,ceil((n)/2)+1);
            else
                v_1 = linspace(a,0,ceil((n)/2));
            end
            v_2 = linspace(0,b,ceil((n)/2));
            v_2 = v_2(2:end); %remove first zero
            v = [v_1,v_2];
        end
    end
    
end


