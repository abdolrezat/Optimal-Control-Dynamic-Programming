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
                this.v_min = -10;
                this.v_max = +10;
                this.n_mesh_v = 10;
                
                this.x_min = -15;
                this.x_max = 15;
                this.n_mesh_x = 10;
                
                %att
                this.w_min = -deg2rad(50);
                this.w_max = -deg2rad(-50);
                this.n_mesh_w = 10;
                
                this.theta1_min = -30; %angles of rotation about y-axis (pitch)
                this.theta1_max = 30;
                this.theta2_min = -20; %angles of rotation about z-axis (yaw)
                this.theta2_max = 20;
                this.theta3_min = -35; %x-axis rotation
                this.theta3_max = 35;
                this.n_mesh_t = 15;
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
                
                this.Qt1 = 6;
                this.Qt2 = 6;
                this.Qt3 = 6;
                this.Qw1 = 6;
                this.Qw2 = 6;
                this.Qw3 = 6;
                
                this.R1 = 4;
                this.R2 = 4;
                this.R3 = 4;
                
                this.T_final = 30;
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
            J_current_1 = J_current_reshaped(obj, s_x1,s_v1,s_t1,s_w1,...
                obj.F_Thr0, obj.F_Thr1, obj.F_Thr6, obj.F_Thr7,...
                obj.Qx1,obj.Qv1,obj.Qt1,obj.Qw1,obj.R1);
            F1 = griddedInterpolant({s_x1,s_v1,s_t1,s_w1},...
                zeros(obj.n_mesh_x,obj.n_mesh_v,obj.n_mesh_t,obj.n_mesh_w,'single'),'linear');
            
            J_current_2 = J_current_reshaped(obj, s_x2,s_v2,s_t2,s_w2,...
                obj.F_Thr2, obj.F_Thr3, obj.F_Thr8, obj.F_Thr9,...
                obj.Qx2,obj.Qv2,obj.Qt2,obj.Qw2,obj.R2);
            F2 = griddedInterpolant({s_x2,s_v2,s_t2,s_w2},...
                zeros(obj.n_mesh_x,obj.n_mesh_v,obj.n_mesh_t,obj.n_mesh_w,'single'),'linear');
            
            J_current_3 = J_current_reshaped(obj, s_x3,s_v3,s_t3,s_w3,...
                obj.F_Thr4, obj.F_Thr5, obj.F_Thr10, obj.F_Thr11,...
                obj.Qx3,obj.Qv3,obj.Qt3,obj.Qw3,obj.R3);
            F3 = griddedInterpolant({s_x3,s_v3,s_t3,s_w3},...
                zeros(obj.n_mesh_x,obj.n_mesh_v,obj.n_mesh_t,obj.n_mesh_w,'single'),'linear');
            
            %calculating next stage states
            %
            [x1_next,v1_next,t1_next,w1_next] = next_stage_states_simplified(obj, s_x1,s_v1,s_t1,s_w1,...
                obj.F_Thr0, obj.F_Thr1, obj.F_Thr6, obj.F_Thr7, obj.J1);
            [x2_next,v2_next,t2_next,w2_next] = next_stage_states_simplified(obj, s_x2,s_v2,s_t2,s_w2,...
                obj.F_Thr2, obj.F_Thr3, obj.F_Thr8, obj.F_Thr9, obj.J2);
            [x3_next,v3_next,t3_next,w3_next] = next_stage_states_simplified(obj, s_x3,s_v3,s_t3,s_w3,...
                obj.F_Thr4, obj.F_Thr5, obj.F_Thr10, obj.F_Thr11, obj.J3);
            
            %beginning (reverse) stage calculations
            whandle = waitbar(0,'Calculation in Progress...');
            for k_s = obj.N_stage-1:-1:1
                tic
                [val, obj.Opt_F_Thr7] = min( J_current_1 + F1(x1_next,v1_next,t1_next,w1_next), [], 8);
                [val, obj.Opt_F_Thr6] = min(  val, [], 7);
                [val, obj.Opt_F_Thr1] = min(  val, [], 6);
                [F1.Values, obj.Opt_F_Thr0] = min( val, [], 5);
                
                [val, obj.Opt_F_Thr9] = min( J_current_2 + F2(x2_next,v2_next,t2_next,w2_next), [], 8);
                [val, obj.Opt_F_Thr8] = min(  val, [], 7);
                [val, obj.Opt_F_Thr3] = min(  val, [], 6);
                [F2.Values, obj.Opt_F_Thr2] = min( val, [], 5);
                
                [val, obj.Opt_F_Thr11] = min( J_current_3 + F3(x3_next,v3_next,t3_next,w3_next), [], 8);
                [val, obj.Opt_F_Thr10] = min(  val, [], 7);
                [val, obj.Opt_F_Thr5] = min(  val, [], 6);
                [F3.Values, obj.Opt_F_Thr4] = min( val, [], 5);
                
                waitbar( 1 - k_s/obj.N_stage, whandle);
                fprintf('step %d - %f seconds\n', k_s, toc)
            end
            keyboard
            
            %get U* Optimal idx
            obj.Opt_F_Thr7 = obj.Opt_F_Thr7(obj.Opt_F_Thr6(obj.Opt_F_Thr1(obj.Opt_F_Thr0)));
            obj.Opt_F_Thr6 = obj.Opt_F_Thr6(obj.Opt_F_Thr1(obj.Opt_F_Thr0));
            obj.Opt_F_Thr1 = obj.Opt_F_Thr1(obj.Opt_F_Thr0);
            obj.Opt_F_Thr0 = obj.Opt_F_Thr0;
            
            obj.Opt_F_Thr9 = obj.Opt_F_Thr9(obj.Opt_F_Thr8(obj.Opt_F_Thr3(obj.Opt_F_Thr2)));
            obj.Opt_F_Thr8 = obj.Opt_F_Thr8(obj.Opt_F_Thr3(obj.Opt_F_Thr2));
            obj.Opt_F_Thr3 = obj.Opt_F_Thr3(obj.Opt_F_Thr2);
            obj.Opt_F_Thr2 = obj.Opt_F_Thr2;
            
            obj.Opt_F_Thr11 = obj.Opt_F_Thr11(obj.Opt_F_Thr10(obj.Opt_F_Thr5(obj.Opt_F_Thr4)));
            obj.Opt_F_Thr10 = obj.Opt_F_Thr10(obj.Opt_F_Thr5(obj.Opt_F_Thr4));
            obj.Opt_F_Thr5 = obj.Opt_F_Thr5(obj.Opt_F_Thr4);
            obj.Opt_F_Thr4 = obj.Opt_F_Thr4;
            
            %set U* Optimal values
            obj.Opt_F_Thr0 = griddedInterpolant({s_x1,s_v1,s_t1,s_w1},...
                obj.F_Thr0(obj.Opt_F_Thr0),...
                'nearest');
            obj.Opt_F_Thr1 = griddedInterpolant({s_x1,s_v1,s_t1,s_w1},...
                obj.F_Thr1(obj.Opt_F_Thr1),...
                'nearest');
            obj.Opt_F_Thr6 = griddedInterpolant({s_x1,s_v1,s_t1,s_w1},...
                obj.F_Thr6(obj.Opt_F_Thr6),...
                'nearest');
            obj.Opt_F_Thr7 = griddedInterpolant({s_x1,s_v1,s_t1,s_w1},...
                obj.F_Thr7(obj.Opt_F_Thr7),...
                'nearest');
            
            
            obj.Opt_F_Thr2 = griddedInterpolant({s_x2,s_v2,s_t2,s_w2},...
                obj.F_Thr2(obj.Opt_F_Thr2),...
                'nearest');
            obj.Opt_F_Thr3 = griddedInterpolant({s_x2,s_v2,s_t2,s_w2},...
                obj.F_Thr3(obj.Opt_F_Thr3),...
                'nearest');
            obj.Opt_F_Thr8 = griddedInterpolant({s_x2,s_v2,s_t2,s_w2},...
                obj.F_Thr8(obj.Opt_F_Thr8),...
                'nearest');
            obj.Opt_F_Thr9 = griddedInterpolant({s_x2,s_v2,s_t2,s_w2},...
                obj.F_Thr9(obj.Opt_F_Thr9),...
                'nearest');
            
            obj.Opt_F_Thr4 = griddedInterpolant({s_x3,s_v3,s_t3,s_w3},...
                obj.F_Thr4(obj.Opt_F_Thr4),...
                'nearest');
            obj.Opt_F_Thr5 = griddedInterpolant({s_x3,s_v3,s_t3,s_w3},...
                obj.F_Thr5(obj.Opt_F_Thr5),...
                'nearest');
            obj.Opt_F_Thr10 = griddedInterpolant({s_x3,s_v3,s_t3,s_w3},...
                obj.F_Thr10(obj.Opt_F_Thr10),...
                'nearest');
            obj.Opt_F_Thr11 = griddedInterpolant({s_x3,s_v3,s_t3,s_w3},...
                obj.F_Thr11(obj.Opt_F_Thr11),...
                'nearest');
            
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
            f2 = reshape(f2,[1 1 1 1 1 length(f2)]  );
            f6 = reshape(f6,[1 1 1 1 1 1 length(f6)]  );
            f7 = reshape(f7,[1 1 1 1 1 1 1 length(f7)]  );
            
            % ODE solve
            x_next = RK4_x(obj, X, V);
            v_next = RK4_v(obj, V, f1,f2,f6,f7);
            t_next = RK4_t(obj, T, W);
            w_next = RK4_w(obj, W, f1,f2,f6,f7 , J);
            

            %repmat each matrix to full size, as required for F inputs
            x_next = repmat(x_next,[1 1 obj.n_mesh_t obj.n_mesh_w length(f1) length(f2) length(f6) length(f7)]);
            v_next = repmat(v_next,[obj.n_mesh_x 1 obj.n_mesh_t obj.n_mesh_w 1 1 1 1]);
            t_next = repmat(t_next,[obj.n_mesh_x obj.n_mesh_v 1 1 length(f1) length(f2) length(f6) length(f7)]);
            w_next = repmat(w_next,[obj.n_mesh_x obj.n_mesh_v obj.n_mesh_t 1 1 1 1 1]);    
            end
        
        function X2 = RK4_x(obj, X1, V)
            % Runge-Kutta - 4th order
            % h = dt;
            k1 = xdynamics(obj, V);
            k2 = xdynamics(obj,(V + k1*obj.h/2));
            k3 = xdynamics(obj,(V + k2*obj.h/2));
            k4 = xdynamics(obj,(V + k3*obj.h));
            
            X2 = X1 + obj.h*(k1 + 2*k2 + 2*k3 + k4)/6;
            
        end
        
        function x_dot = xdynamics(~,v)
            x_dot = v;
        end
        
        function V2 = RK4_v(obj, V1, f1,f2,f6,f7) % does not need RK4, ki's are equal
            % Runge-Kutta - 4th order
            % h = dt;
            k1 = vdynamics(obj, V1 , f1,f2,f6,f7);
            k2 = vdynamics(obj,(V1 + k1*obj.h/2), f1,f2,f6,f7);
            k3 = vdynamics(obj,(V1 + k2*obj.h/2), f1,f2,f6,f7);
            k4 = vdynamics(obj,(V1 + k3*obj.h), f1,f2,f6,f7);
            
            V2 = V1 + obj.h*(k1 + 2*k2 + 2*k3 + k4)/6;
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
            k2 = tdynamics(obj,(W1 + k1*obj.h/2));
            k3 = tdynamics(obj,(W1 + k2*obj.h/2));
            k4 = tdynamics(obj,(W1 + k3*obj.h));
            
            T2 = T1 + obj.h*(k1 + 2*k2 + 2*k3 + k4)/6;
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
            k2 = wdynamics(obj,(w + k1*obj.h/2), f1,f2,f6,f7, J);
            k3 = wdynamics(obj,(w + k2*obj.h/2), f1,f2,f6,f7, J);
            k4 = wdynamics(obj,(w + k3*obj.h), f1,f2,f6,f7, J);
            
            W2 = w + obj.h*(k1 + 2*k2 + 2*k3 + k4)/6;
        end
        
        
        function w_dot = wdynamics(obj,~, f1,f2,f6,f7, J)
            w_dot = (f1*obj.T_dist +f2*(-obj.T_dist) + ...
                f6*obj.T_dist +f7*(-obj.T_dist) )/J;
        end
        
        function [f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11] = get_thruster_on_off_optimal(~,x,v,t,w)
            t_x = t(1); %rotation about x_axis
            t_y = t(2);
            t_z = t(3);
            
            w_x = w(1); %rotational speed about x_axis
            w_y = w(2);
            w_z = w(3);
            
            f0 = obj.Opt_F_Thr0(x(1),v(1), t_y, w_y );
            f1 = obj.Opt_F_Thr1(x(1),v(1), t_y, w_y );
            f6 = obj.Opt_F_Thr6(x(1),v(1), t_y, w_y );
            f7 = obj.Opt_F_Thr7(x(1),v(1), t_y, w_y );
            
            f2 = obj.Opt_F_Thr2(x(2),v(2), t_z, w_z );
            f3 = obj.Opt_F_Thr3(x(2),v(2), t_z, w_z );
            f8 = obj.Opt_F_Thr8(x(2),v(2), t_z, w_z );
            f9 = obj.Opt_F_Thr9(x(2),v(2), t_z, w_z );
            
            f4 = obj.Opt_F_Thr4(x(3),v(3), t_x, w_x );
            f5 = obj.Opt_F_Thr5(x(3),v(3), t_x, w_x );
            f10 = obj.Opt_F_Thr10(x(3),v(3), t_x, w_x );
            f11 = obj.Opt_F_Thr11(x(3),v(3), t_x, w_x );
            
        end
        
        
        function get_optimal_path(obj)
            X0 = obj.defaultX0;
            tspan = 0:obj.h:obj.T_final;
            X_ode45 = zeros(obj.N_stage, 13);
            X_ode45(1,:) = X0;
            
            for k_stage=1:obj.N_stage-1
                %determine F_Opt each Thruster
                x = X_ode45(k_stage,1:3);
                t = [2*asin(X_stage(6));... %angle about x-axis
                    2*asin(X_stage(5));... %y-axis
                    2*asin(X_stage(4))]; %z-axis
                [f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11] = get_thruster_on_off_optimal(obj,x,v,t,w)
                %
                [~,X_temp] = ode45(@ode_eq,[tspan(k_stage), tspan(k_stage+1)], X_ode45(k_stage,:));
                X_ode45(k_stage+1,:) = X_temp(end,:);
            end
            
            %function declarations
            function x_dot = ode_eq(~,X1)
                x_dot = system_dynamics(X1);
                x_dot = x_dot';
                function X_dot = system_dynamics(X, U)
                    x1 = X(1);
                    x2 = X(2);
                    x3 = X(3);
                    
                    q1 = X(4);
                    q2 = X(5);
                    q3 = X(6);
                    q4 = X(7);
                    
                    %InertiaM is a complete matrix
                    w = [x1;x2;x3];
                    w_dot = obj.InertiaM\(U - cross(w, obj.InertiaM*w));
                    
                    X_dot(1) = w_dot(1);
                    X_dot(2) = w_dot(2);
                    X_dot(3) = w_dot(3);
                    
                    X_dot(4) = 0.5*(x3.*q2 -x2.*q3 +x1.*q4);
                    X_dot(5) = 0.5*(-x3.*q1 +x1.*q3 +x2.*q4);
                    X_dot(6) = 0.5*(x2.*q1 -x1.*q2 +x3.*q4);
                    X_dot(7) = 0.5*(-x1.*q1 -x2.*q2 -x3.*q3);
                end
            end
        end
        
        function get_optimal_path_pos(obj)
            %%
            global mu
            mu = 398600;
            
            %   Prescribed initial state vector of chaser B in the co-moving frame:
            dr0 = [-1  0  0];
            dv0 = [ 0 0 0];
            y0 = [dr0 dv0]';
            
            tf  = obj.T_final;
            %...End input data
            
            %...Calculate the target's initial state vector using Algorithm 4.5:
            [R0,V0] = get_target_R0V0(obj);
            
            %%
            N = ceil(tf/obj.h);
            tspan = 0:obj.h:tf;
            X_ode45 = zeros(6, N);
            F_Opt_history = zeros(3, N);
            X_ode45(:,1) = y0;
            tic
            for k_stage=1:N-1
                %determine U from X_stage
                X_stage = X_ode45(:,k_stage);
                a_x = obj.U1_Opt(X_stage(1), X_stage(4)); %x1,v1
                a_y = obj.U2_Opt(X_stage(2), X_stage(5)); %x2,v2
                a_z = obj.U3_Opt(X_stage(3), X_stage(6)); %x3,v3
                F_Opt_history(:,k_stage) = [a_x;a_y;a_z];
                % variables used in nested differential equation function
                % moved to inside @rates
                
                %
                [~,X_temp] = rkf45(@rates,[tspan(k_stage), tspan(k_stage+1)], X_ode45(:,k_stage));
                X_ode45(:,k_stage+1) = X_temp(end,:)';
            end
            toc
            
            T_ode45 = tspan(1:end-1)';
            
            %plot states x
            figure
            hold on
            grid on
            for i=1:3
                plot(T_ode45,X_ode45(i,:))
            end
            legend('x1','x2','x3')
            
            %plot states v
            figure
            hold on
            grid on
            for i=4:6
                plot(T_ode45,X_ode45(i,:))
            end
            legend('v1','v2','v3')
            
            %plot controls
            figure
            hold on
            grid on
            for i=1:3
                plot(T_ode45,F_Opt_history(i,:))
            end
            legend('u1','u2','u3')
            
            
            
            
            function dydt = rates(t,y)
                % ~~~~~~~~~~~~~~~~~~~~~~~~
                %{
  This function computes the components of f(t,y) in Equation 7.36.
  
  t             - time
  f             - column vector containing the relative position and
                  velocity vectors of B at time t
  R, V          - updated state vector of A at time t
  X, Y, Z       - components of R
  VX, VY, VZ    - components of V
  R_            - magnitude of R
  RdotV         - dot product of R and V
  h             - magnitude of the specific angular momentum of A

  dx , dy , dz  - components of the relative position vector of B
  dvx, dvy, dvz - components of the relative velocity vector of B
  dax, day, daz - components of the relative acceleration vector of B
  dydt          - column vector containing the relative velocity
                  and acceleration components of B at time t

  User M-function required: rv_from_r0v0
                %}
                % ------------------------
                %...Update the state vector of the target orbit using Algorithm 3.4:
                
                %                 X  = R(1); Y  = R(2); Z  = R(3);
                %                 VX = V(1); VY = V(2); VZ = V(3);
                %
                [R,V] = update_RV_target(obj, R0, V0, t);
                norm_R = (R*R')^.5; %norm R
                RdotV = sum(R.*V); %dot product
                crossRV = [R(2).*V(3)-R(3).*V(2); % cross product of R and V
                    R(3).*V(1)-R(1).*V(3);
                    R(1).*V(2)-R(2).*V(1)];
                H  = (crossRV'*crossRV)^.5 ; %norm(crossRV);
                
                dx  = y(1);
                dy  = y(2);
                dz  = y(3);
                
                dvx = y(4);
                dvy = y(5);
                dvz = y(6);
                
                dax   =  (2*mu/norm_R^3 + H^2/norm_R^4)*dx - 2*RdotV/norm_R^4*H*dy + 2*H/norm_R^2*dvy + a_x;
                day   =   -(mu/norm_R^3 - H^2/norm_R^4)*dy + 2*RdotV/norm_R^4*H*dx - 2*H/norm_R^2*dvx + a_y;
                daz   = -mu/norm_R^3*dz + a_z;
                
                dydt  = [dvx dvy dvz dax day daz]';
            end %rates
            
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
            f2 = reshape(f2,[1 1 1 1 1 length(f2)]  );
            f3 = reshape(f3,[1 1 1 1 1 1 length(f3)]  );
            f4 = reshape(f4,[1 1 1 1 1 1 1 length(f4)]  );
            
            J_current_M = (Qx * x.^2 + Qv * v.^2 + Qw * w.^2 + Qt * t.^2 +...
                ( R* f1.^2 +  R* f2.^2 +  R* f3.^2 +  R* f4.^2)  );
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


