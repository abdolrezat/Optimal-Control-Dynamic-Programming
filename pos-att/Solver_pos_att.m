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
        
        T_dist %Thruster placement distance
        M_Thr0 %Thruster 0 - Moment
        M_Thr1
        M_Thr2
        M_Thr3
        M_Thr4
        M_Thr5
        M_Thr6
        M_Thr7
        M_Thr8
        M_Thr9
        M_Thr10
        M_Thr11
        %
        %Optimum values for each Thruster
        Opt_F_Thr0 single
        Opt_F_Thr1 single
        Opt_F_Thr2 single
        Opt_F_Thr3 single
        Opt_F_Thr4 single
        Opt_F_Thr5 single
        Opt_F_Thr6 single
        Opt_F_Thr7 single
        Opt_F_Thr8 single
        Opt_F_Thr9 single
        Opt_F_Thr10 single
        Opt_F_Thr11 single

        
    end
    
    methods
        function this = Solver_position()
            if nargin < 1
                
                %pos
                this.v_min = -10;
                this.v_max = +10;
                this.n_mesh_v = 300;
                
                this.x_min = -15;
                this.x_max = 15;
                this.n_mesh_x = 100;
                
                %att
                this.w_min = -deg2rad(50);
                this.w_max = -deg2rad(-50);
                this.n_mesh_w = 100;
                
                this.t1_min = -30; %angles of rotation about y-axis (pitch)
                this.t1_max = 30;
                this.t2_min = -20; %angles of rotation about z-axis (yaw)
                this.t2_max = 20;
                this.t3_min = -35; %x-axis rotation
                this.t3_max = 35;
                this.n_mesh_t = 300;
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
            
            this.defaultX0 = [0;0;0;0;0;0]; % [x;v]
            
            %Thruster Forces, Moments
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
            
            %y-direction
            this.M_Thr0 =  this.F_Thr0*this.T_dist;
            this.M_Thr1 = -this.F_Thr1*this.T_dist;
            
            this.M_Thr6 =  this.F_Thr6*this.T_dist;
            this.M_Thr7 = -this.F_Thr7*this.T_dist;
            
            %z-direction
            this.M_Thr2 =  this.F_Thr2*this.T_dist;
            this.M_Thr3 = -this.F_Thr3*this.T_dist;
            
            this.M_Thr8 =  this.F_Thr8*this.T_dist;
            this.M_Thr9 = -this.F_Thr9*this.T_dist;
            
            %x-direction
            this.M_Thr4 =  this.F_Thr4*this.T_dist;
            this.M_Thr5 = -this.F_Thr5*this.T_dist;
            
            this.M_Thr10 =  this.F_Thr10*this.T_dist;
            this.M_Thr11 = -this.F_Thr11*this.T_dist;
            
            
            
        end
        
        function simplified_run(obj)
            
            %% mesh generation
            s_x1 = linspace(obj.x_min, obj.x_max, obj.n_mesh_x);
            s_x2 = linspace(obj.x_min, obj.x_max, obj.n_mesh_x);
            s_x3 = linspace(obj.x_min, obj.x_max, obj.n_mesh_x);
            
            s_v1 = linspace(obj.v_min, obj.v_max, obj.n_mesh_v);
            s_v2 = linspace(obj.v_min, obj.v_max, obj.n_mesh_v);
            s_v3 = linspace(obj.v_min, obj.v_max, obj.n_mesh_v);
            
            s_t1 = linspace(deg2rad(obj.theta1_min), deg2rad(obj.theta1_max), obj.n_mesh_t);
            s_t2 = linspace(deg2rad(obj.theta2_min), deg2rad(obj.theta2_max), obj.n_mesh_t);
            s_t3 = linspace(deg2rad(obj.theta3_min), deg2rad(obj.theta3_max), obj.n_mesh_t);
            
            s_w1 = linspace(obj.w_min, obj.w_max, obj.n_mesh_w);
            s_w2 = linspace(obj.w_min, obj.w_max, obj.n_mesh_w);
            s_w3 = linspace(obj.w_min, obj.w_max, obj.n_mesh_w);
            
            
            %% initialization
            %states, controls
            %% issue: change T,W with the appropriate theta/w that are affected by the 4 Thrusters
            [X1,V1,T1,W1,Fv0,Fv1,Fv6,Fv7] = ndgrid(s_x1,s_v1,s_t1,s_w1,... %% moment in y-direction ***
                this.F_Thr0, this.F_Thr1, this.F_Thr6, this.F_Thr7);
            [X2,V2,T2,W2,Fv2,Fv3,Fv8,Fv9] = ndgrid(s_x2,s_v2,s_t2,s_w2,...
                this.F_Thr2, this.F_Thr3, this.F_Thr8, this.F_Thr9);
            [X3,V3,T3,W3,Fv4,Fv5,Fv10,Fv11] = ndgrid(s_x3,s_v3,s_t3,s_w3,...
                this.F_Thr4, this.F_Thr5, this.F_Thr10, this.F_Thr11);
            
            J_current = @(x,v,t,w,f1,f2,f3,f4,Qx,Qv,Qt,Qw,R)...
                (Qx * x.^2 + Qv * v.^2 + Qw * w.^2 + Qt * t.^2 +...
                R* (f1.^2 + f2.^2 + f3.^2 + f4.^2)  );
            %calculating J fixed
            J_current_1 = J_current(X1,V1,T1,W1,Fv0,Fv1,Fv6,Fv7,obj.Qx1,obj.Qv1,obj.Qt1,obj.Qw1,obj.R1);
            F1 = griddedInterpolant({s_x1,s_v1,s_t1,s_w1},...
                zeros(obj.n_mesh_x,obj.n_mesh_v,obj.n_mesh_t,obj.n_mesh_w,'single'),'linear');
            
            J_current_2 = J_current(X2,V2,T2,W2,Fv2,Fv3,Fv8,Fv9,obj.Qx2,obj.Qv2,obj.Qt2,obj.Qw2,obj.R2);
            F2 = griddedInterpolant({s_x2,s_v2,s_t2,s_w2},...
                zeros(obj.n_mesh_x,obj.n_mesh_v,obj.n_mesh_t,obj.n_mesh_w,'single'),'linear');
            
            J_current_3 = J_current(X3,V3,T3,W3,Fv4,Fv5,Fv10,Fv11,obj.Qx3,obj.Qv3,obj.Qt3,obj.Qw3,obj.R3);
            F3 = griddedInterpolant({s_x3,s_v3,s_t3,s_w3},...
                zeros(obj.n_mesh_x,obj.n_mesh_v,obj.n_mesh_t,obj.n_mesh_w,'single'),'linear');
            
            %calculating next stage states
            %
            [x1_next,v1_next,t1_next,w1_next] = next_stage_states(obj, X1, V1, T1, W1, Fv0,Fv1,Fv6,Fv7  , J);
            [x2_next,v2_next,t2_next,w2_next] = next_stage_states(obj, X2, V2, T2, W2, Fv2,Fv3,Fv8,Fv9  , J);
            [x3_next,v3_next,t3_next,w3_next] = next_stage_states(obj, X3, V3, T3, W3, Fv4,Fv5,Fv10,Fv11, J);
            
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
        
        function [x_next,v_next,t_next,w_next] = next_stage_states(obj, X, V, T, W, f1,f2,f6,f7 ,J)
            x_next = RK4_x(obj, X, V);
            v_next = RK4_v(obj, V, f1,f2,f6,f7);
            t_next = RK4_t(obj, T, W);
            w_next = RK4_w(obj, W, f1,f2,f6,f7 , J);
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
        
        
        function get_optimal_path_simplified()
            
        end
    end
    
end


