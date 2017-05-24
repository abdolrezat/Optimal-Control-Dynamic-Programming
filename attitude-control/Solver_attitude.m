classdef Solver_attitude < handle
    %SOLVER_ATTITUDE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        N % number of stages
        InertiaM % Moment of Inertia Matrix
        J1 % element 1 of Inertia Matrix
        J2 % ...
        J3 % ..
        size_state_mat
        dim_U1 % dimension number of U, after the 7 state variables
        dim_U2
        dim_U3
        U_vector % values that are applied as control inputs
        % big data
        %J_next_states_opt single  % optimum values of J will be stored in this matrix and recalled each stage
        J_current_state_fix single  % values of cost to reach next stage, always fixed in Time Invariant system
        F % gridded Interpolant holding values of "J_next_states_opt"
        %obsolete, use X1V,etc... instead
        %         X1 single % (s+c)-dimensional grid of state 1
        %         X2 single
        %         X3 single
        %         X4 single
        %         X5 single
        %         X6 single
        %         X7 single
        %         U1 single  % (s+c)-dimensional grid of control input 1
        %         U2 single
        %         U3 single
        X1_next single  % (s+c)-dim grid of next state-1
        X2_next single
        X3_next single
        %         X4_next single
        X4_next single
        X5_next single
        X6_next single
        
        U1_Opt   % (s-dim x k) grid of optimal values U* for every stage and states
        U2_Opt
        U3_Opt
        %-- ranges for grid generation
        w_min
        w_max
        n_mesh_w
        n_mesh_t % for simplified dyn
        yaw_min
        yaw_max
        pitch_min
        pitch_max
        roll_min
        roll_max
        n_mesh_q
        sr_1 % state 1 range vector for grid and interpolant generation
        sr_2 % "
        sr_3 % "
        %         sr_4 % "
        s_yaw % state 5
        s_pitch % " 6
        s_roll % " 7
        cr_1 % control input 1 range vector for grid generation
        cr_2
        cr_3
        
        T_final
        h %time step for discrete time system
        N_stage
        
        Q1
        Q2
        Q3
        Q4
        Q5
        Q6
        R1
        R2
        R3
        
        Qt1 %for simplified dyn
        Qt2
        Qt3
        
        X1V single
        X2V single
        X3V single
        cang_x4 single
        sang_x4 single
        cang_x5 single
        sang_x5 single
        cang_x6 single
        sang_x6 single
        cang_x7 single
        sang_x7 single
        U1V single
        U2V single
        U3V single
        
        defaultX0 % default value for initial states if not explicitly specified
    end
    
    methods
        %Constructor
        function this = Solver_attitude()
            % states are [w1; w2; w3; q1; q2; q3; q4]
            if nargin < 1
                w_min = -deg2rad(0.9);
                w_max = -deg2rad(-0.9);
                n_mesh_w = 100;
                this.yaw_min = -30; %angles
                this.yaw_max = 30;
                this.pitch_min = -20;
                this.pitch_max = 20;
                this.roll_min = -35;
                this.roll_max = 35;
                this.n_mesh_q = 10;
                this.n_mesh_t = 300; %for simplified dyn
                
                this.InertiaM = [2  0   0;...
                                 0 2.5  0;...
                                 0  0  3];
                    
                
                this.Q1 = 1;
                this.Q2 = 1;
                this.Q3 = 1;
                
                this.Q4 = 6; %quaternions, thetas
                this.Q5 = 6;
                this.Q6 = 6;
                this.R1 = 4;
                this.R2 = 4;
                this.R3 = 4;
                
                this.Qt1 = this.Q4;
                this.Qt2 = this.Q5;
                this.Qt3 = this.Q6;
                
                this.T_final = 50;
                this.h = 0.005;
            end
            this.w_min = w_min;
            this.w_max = w_max;
            this.n_mesh_w = n_mesh_w;
            this.N_stage = this.T_final/this.h;
            
            if(~isinteger( this.N_stage))
                this.N_stage = ceil(this.N_stage);
                this.T_final = this.h*this.N_stage;
                warning('T_final is not a factor of h (dt), increasing T_final to %.2f\n',this.T_final)
            end
            
            q0 = [0;0;0.0174524064372835;0.999847695156391]; %angle2quat(0,0,roll = deg2rad(2))
            this.defaultX0 = [0;0;0;...
                ...  %quaternions equal to quat(deg2rad(-10), deg2rad(20), deg2rad(-15))
                q0];%0.999660006156261;0.00841930262082080;0.0176013597667272;0.0172968080698774];
            
            
            this.J1 = this.InertiaM(1);
            this.J2 = this.InertiaM(5);
            this.J3 = this.InertiaM(9);
            
            this.dim_U1 = 7;
            this.dim_U2 = 8;
            this.dim_U3 = 9;
            this.U_vector = [-0.01 0 0.01];
            length_U = length(this.U_vector);
            %Preallocation and mesh generation
            this.sr_1 = linspace(w_min, w_max, n_mesh_w);
            this.sr_2 = linspace(w_min, w_max, n_mesh_w);
            this.sr_3 = linspace(w_min, w_max, n_mesh_w);
            
            this.s_yaw = linspace(deg2rad(this.yaw_min), deg2rad(this.yaw_max), this.n_mesh_q);
            this.s_pitch = linspace(deg2rad(this.pitch_min), deg2rad(this.pitch_max), this.n_mesh_q);
            this.s_roll = linspace(deg2rad(this.roll_min), deg2rad(this.roll_max), this.n_mesh_q);
            
            this.size_state_mat = [n_mesh_w,n_mesh_w,n_mesh_w,...
                this.n_mesh_q,this.n_mesh_q,this.n_mesh_q];
            
            this.U1_Opt = zeros([this.size_state_mat],'uint8');
            this.U2_Opt = zeros([this.size_state_mat,length_U],'uint8');
            this.U3_Opt = zeros([this.size_state_mat,length_U,length_U],'uint8');
            
        end
        
        %Calculation of optimal matrices
        function simplified_run(obj)
            
            %% mesh generation
            s_w1 = linspace(obj.w_min, obj.w_max, obj.n_mesh_w);
            s_w2 = linspace(obj.w_min, obj.w_max, obj.n_mesh_w);
            s_w3 = linspace(obj.w_min, obj.w_max, obj.n_mesh_w);
            
            s_t1 = linspace(deg2rad(obj.yaw_min), deg2rad(obj.yaw_max), obj.n_mesh_t);
            s_t2 = linspace(deg2rad(obj.pitch_min), deg2rad(obj.pitch_max), obj.n_mesh_t);
            s_t3 = linspace(deg2rad(obj.roll_min), deg2rad(obj.roll_max), obj.n_mesh_t);
            
            %% initialization
            %states, controls
            [X1,X4,U1] = ndgrid(s_w1,s_t1,obj.U_vector);
            [X2,X5,U2] = ndgrid(s_w2,s_t2,obj.U_vector);
            [X3,X6,U3] = ndgrid(s_w3,s_t3,obj.U_vector);
            
            %u_opt
            size_Umat = [obj.n_mesh_w, obj.n_mesh_t];
            obj.U1_Opt = zeros(size_Umat,'single');
            obj.U2_Opt = zeros(size_Umat,'single');
            obj.U3_Opt = zeros(size_Umat,'single');
            
            %J matrices, current and optimal next stage
            J_current = @(w,theta,U,Qw,Qt,R)(Qw * w.^2 + Qt * theta.^2 + R*U.^2);
            J_current_1 = J_current(X1,X4,U1,obj.Q1,obj.Qt1,obj.R1);
            F1 = griddedInterpolant({s_w1,s_t1},zeros(obj.n_mesh_w,obj.n_mesh_t),'linear');
            
            J_current_2 = J_current(X2,X5,U2,obj.Q2,obj.Qt2,obj.R2);
            F2 = griddedInterpolant({s_w2,s_t2},zeros(obj.n_mesh_w,obj.n_mesh_t),'linear');
            
            J_current_3 = J_current(X3,X6,U3,obj.Q3,obj.Qt3,obj.R3);
            F3 = griddedInterpolant({s_w3,s_t3},zeros(obj.n_mesh_w,obj.n_mesh_t),'linear');
            %
            %
            [w1_next,t1_next] = next_stage_states_simplified(obj,X1, X4, U1, obj.J1, obj.h);
            [w2_next,t2_next] = next_stage_states_simplified(obj,X2, X5, U2, obj.J2, obj.h);
            [w3_next,t3_next] = next_stage_states_simplified(obj,X3, X6, U3, obj.J3, obj.h);
            % keyboard
            whandle = waitbar(0,'Calculation in Progress...');
            for k_s = obj.N_stage-1:-1:1
                tic
                %     C = F1(w1_next,t1_next);
                [F1.Values, U1_idx] = min( J_current_1 + F1(w1_next,t1_next), [], 3);
                [F2.Values, U2_idx] = min( J_current_2 + F2(w2_next,t2_next), [], 3);
                [F3.Values, U3_idx] = min( J_current_3 + F3(w3_next,t3_next), [], 3);
                
                
                %     QQ = U1_idx;
                waitbar( 1 - k_s/obj.N_stage, whandle);
                fprintf('step %d - %f seconds\n', k_s, toc)
            end
            %set U* Optimal
            obj.U1_Opt = griddedInterpolant({s_w1,s_t1}, obj.U_vector(U1_idx),'nearest');
            obj.U2_Opt = griddedInterpolant({s_w2,s_t2}, obj.U_vector(U2_idx),'nearest');
            obj.U3_Opt = griddedInterpolant({s_w3,s_t3}, obj.U_vector(U3_idx),'nearest');
            close(whandle)
            fprintf('stage calculation complete... cleaning up...\n')
            % clear J_current_1 J_current_2 J_current_3 F1 F2 F3 X1 X2 X3 X4 X5 X6 U1 U2 U3 ...
            %     w1_next w2_next w3_next t1_next t2_next t3_next
            fprintf('...Done!\n')
            
            
        end
        
        function run(obj)
            obj.reshape_states();
            % calculate cost to reach next stage
            %obj.J_current_state_fix = g_D(obj);
            fprintf('calculating fixed cost matrix...\n')
            calculate_J_current_state_fix_shaped(obj);
            % Calculate and store J*NN = h(xi(N)) for all x(N)
            obj.F = griddedInterpolant(...
                {obj.sr_1, obj.sr_2, obj.sr_3, obj.s_yaw, obj.s_pitch, obj.s_roll}, ...
                zeros(obj.size_state_mat,'single'),'linear');
            
            %
            fprintf('calculating next stage states...\n')
            calculate_states_next(obj);
            
            %
            %% obj.stage_J_star = H*X;
            %%
            fprintf('beginning stage calculation...\n')
            for k_s=obj.N_stage-1:-1:1
                tic
                calculate_J_U_opt_state_M(obj, k_s);
                fprintf('step %d - %f seconds\n', k_s, toc)
                %stop criteria (e.g. dFi(X) < tol, etc...) can be added here
                %                 id1 = this.sr_1 ==
                %
            end
            %% final Ui_Opts
            % do NOT reverse the order of assignments
            obj.U3_Opt = obj.U3_Opt(obj.U2_Opt(obj.U1_Opt));
            obj.U2_Opt = obj.U2_Opt(obj.U1_Opt);
            obj.U1_Opt = obj.U1_Opt;
            %report
            fprintf('stage calculation complete... cleaning up\n')
            clear_up(obj)
            obj.U1_Opt = single(obj.U_vector(obj.U1_Opt));
            obj.U2_Opt = single(obj.U_vector(obj.U2_Opt));
            obj.U3_Opt = single(obj.U_vector(obj.U3_Opt));
            fprintf('...Done!\n')
        end
        
        function [X1,X2,X3,X4,X5,X6,X7] = a_D_M(obj, X1,X2,X3,X4,X5,X6,X7,U1,U2,U3)
            
            
        end
        
        function calculate_J_current_state_fix(obj)
            obj.J_current_state_fix = obj.Q1*obj.X1.^2 + obj.Q2*obj.X2.^2 + ...
                obj.Q3*obj.X3.^2 + ...  % obj.Q4*obj.X4.^2 + ...
                obj.Q5*obj.X5.^2 + obj.Q6*obj.X6.^2 + ...
                obj.Q7*obj.X7.^2 + obj.R1*obj.U1.^2 + ...
                obj.R2*obj.U2.^2 + obj.R3*obj.U3.^2;
        end
        
        function calculate_J_current_state_fix_shaped(obj)
            obj.J_current_state_fix = obj.Q1*obj.X1V.^2 + obj.Q2*obj.X2V.^2 + ...
                obj.Q3*obj.X3V.^2 + ...  % obj.Q4*obj.X4.^2 + ...
                obj.Q4*(obj.sang_x4.*obj.cang_x5.*obj.cang_x6 - obj.cang_x4.*obj.sang_x5.*obj.sang_x6).^2 + ...
                obj.Q5*(obj.cang_x4.*obj.sang_x5.*obj.cang_x6 + obj.sang_x4.*obj.cang_x5.*obj.sang_x6).^2 + ...
                obj.Q6*(obj.cang_x4.*obj.cang_x5.*obj.sang_x6 - obj.sang_x4.*obj.sang_x5.*obj.cang_x6).^2 + ...
                obj.R1*obj.U1V.^2 + obj.R2*obj.U2V.^2 + obj.R3*obj.U3V.^2;
            %%
            %q's are defined as follows:
            %x4 = q1 as defined in Kirk Control Theory, q4 by MATLAB
            %definition
            % (obj.sang_x4.*obj.cang_x5.*obj.cang_x6 - obj.cang_x4.*obj.sang_x5.*obj.sang_x6)
            
            %x5 = q2 as defined in Kirk Control Theory, q3 by MATLAB
            %definition
            % (obj.cang_x4.*obj.sang_x5.*obj.cang_x6 + obj.sang_x4.*obj.cang_x5.*obj.sang_x6)
            
            %x6 = q3 as defined in Kirk Control Theory, q2 by MATLAB
            %definition
            %%(obj.cang_x4.*obj.cang_x5.*obj.sang_x6 - obj.sang_x4.*obj.sang_x5.*obj.cang_x6)
            
            %x7 = q4 as defined in Kirk Control Theory, q1 by MATLAB
            %definition
            %(obj.cang_x4.*obj.cang_x5.*obj.cang_x6 + obj.sang_x4.*obj.sang_x5.*obj.sang_x6)
            
            %%
            
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
        
        function J = g_D(obj,X1,X2,X3,X5,X6,X7,U1,U2,U3)
            J = obj.Q1*X1.^2 + obj.Q2*X2.^2 + ...
                obj.Q3*X3.^2 + ... %obj.Q4*X4.^2 + ...
                obj.Q5*X5.^2 + obj.Q6*X6.^2 + ...
                obj.Q7*X7.^2 + obj.R1*U1.^2 + ...
                obj.R2*U2.^2 + obj.R3*U3.^2;
        end
        
        function calculate_J_U_opt_state_M(obj)
            
            %% CAUTION: this interpolant is only valid for Xmesh
            %find J final for each state and control (X,U) and add it to next state
            %optimum J*
            %             nq = obj.n_mesh_q;
            %             nw = obj.n_mesh_w;
            %             nu = length(obj.U_vector);
            %             [val, U_ID3] = min( obj.J_current_state_fix  + ...
            %                 F(repmat(obj.X1_next,[1 1 1 nq nq nq 1 nu nu]),... %X1_next size = nw x nw x nw x 1 x 1 x 1 x nu
            %                 repmat(obj.X2_next,[1 1 1 nq nq nq nu 1 nu]),... %X2_next size = nw x nw x nw x 1 x 1 x 1 x 1 x nu
            %                 repmat(obj.X3_next,[1 1 1 nq nq nq nu nu 1]),... %X3_next size = nw x nw x nw x 1 x 1 x 1 x 1 x 1 x nu
            %                 repmat(obj.X5_next,[1 1 1 1 1 1 nu nu nu]),... %X5_next size = nw x nw x nw x nq x nq x nq
            %                 repmat(obj.X6_next,[1 1 1 1 1 1 nu nu nu]),... %X6_next size = nw x nw x nw x nq x nq x nq
            %                 repmat(obj.X7_next,[1 1 1 1 1 1 nu nu nu])), ... %X7_next size = nw x nw x nw x nq x nq x nq
            %                 [], obj.dim_U3);
            [val, obj.U3_Opt] = min( obj.J_current_state_fix  + ...
                obj.F(obj.X1_next,... %X1_next size = nw x nw x nw x 1 x 1 x 1 x nu
                obj.X2_next,... %X2_next size = nw x nw x nw x 1 x 1 x 1 x 1 x nu
                obj.X3_next,... %X3_next size = nw x nw x nw x 1 x 1 x 1 x 1 x 1 x nu
                obj.X4_next,... %X4_next size = nw x nw x nw x nq x nq x nq
                obj.X5_next,... %X5_next size = nw x nw x nw x nq x nq x nq
                obj.X6_next), ... %X6_next size = nw x nw x nw x nq x nq x nq
                [], obj.dim_U3);
            [val, obj.U2_Opt] = min( val, [], obj.dim_U2);
            [obj.F.Values, obj.U1_Opt] = min( val, [], obj.dim_U1);
            
        end
        
        function spacecraft_dynamics_taylor_estimate(obj)
            %returns the derivatives x_dot = f(X,u)
            %J is assumed to be diagonal, J12 = J23 = ... = 0
            %             x4 = (1 - (q1.^2 + ...
            %                 q2.^2 + ...
            %                 q3.^2)).^0.5;
            x7 = (1 - ( (obj.sang_x4.*obj.cang_x5.*obj.cang_x6 - obj.cang_x4.*obj.sang_x5.*obj.sang_x6).^2 + ...
                (obj.cang_x4.*obj.sang_x5.*obj.cang_x6 + obj.sang_x4.*obj.cang_x5.*obj.sang_x6).^2 + ...
                (obj.cang_x4.*obj.cang_x5.*obj.sang_x6 - obj.sang_x4.*obj.sang_x5.*obj.cang_x6).^2) ).^0.5;
            
            obj.X1_next = obj.X1V + obj.h*((obj.J2-obj.J3)/obj.J1*obj.X2V.*obj.X3V   +  obj.U1V/obj.J1);
            obj.X2_next = obj.X2V + obj.h*((obj.J3-obj.J1)/obj.J2*obj.X3V.*obj.X1V   +  obj.U2V/obj.J2);
            obj.X3_next = obj.X3V + obj.h*((obj.J1-obj.J2)/obj.J3*obj.X1V.*obj.X2V   +  obj.U3V/obj.J3);
            %obj.X4_next = obj.X4 + obj.h*(0.5*(-obj.X1.*obj.X7 -obj.X2.*obj.X6 -obj.X3.*obj.X5));
            % simplified:
            %             obj.X4_next = x4 + ...
            %                 obj.h*(0.5*(obj.X3V.* x5 ...
            %                 -obj.X2V.*x6 ...
            %                 +obj.X1V.*x7));
            %
            %             obj.X5_next = x5 + ...
            %                 obj.h*(0.5*(-obj.X3V.*x4 ...
            %                 +obj.X1V.*x6 ...
            %                 +obj.X2V.*x7));
            %
            %             obj.X6_next = x6 + ...
            %                 obj.h*(0.5*(obj.X2V.*x4 ...
            %                 -obj.X1V.*x5 ...
            %                 +obj.X3V.*x7));
            %
            %             %x4_next
            %             x7 = x7 + obj.h*(0.5*(-obj.X1V.*x4 ...
            %                 -obj.X2V.*x5 ...
            %                 -obj.X3V.*x6 ));
            
            %             % expanded(actual):
            obj.X4_next = (obj.sang_x4.*obj.cang_x5.*obj.cang_x6 - obj.cang_x4.*obj.sang_x5.*obj.sang_x6) + ...
                obj.h*(0.5*(obj.X3V.* (obj.cang_x4.*obj.sang_x5.*obj.cang_x6 + obj.sang_x4.*obj.cang_x5.*obj.sang_x6) ...
                -obj.X2V.*(obj.cang_x4.*obj.cang_x5.*obj.sang_x6 - obj.sang_x4.*obj.sang_x5.*obj.cang_x6) ...
                +obj.X1V.*x7));
            
            obj.X5_next = (obj.cang_x4.*obj.sang_x5.*obj.cang_x6 + obj.sang_x4.*obj.cang_x5.*obj.sang_x6) + ...
                obj.h*(0.5*(-obj.X3V.*(obj.sang_x4.*obj.cang_x5.*obj.cang_x6 - obj.cang_x4.*obj.sang_x5.*obj.sang_x6) ...
                +obj.X1V.*(obj.cang_x4.*obj.cang_x5.*obj.sang_x6 - obj.sang_x4.*obj.sang_x5.*obj.cang_x6) ...
                +obj.X2V.*x7));
            
            obj.X6_next = (obj.cang_x4.*obj.cang_x5.*obj.sang_x6 - obj.sang_x4.*obj.sang_x5.*obj.cang_x6) + ...
                obj.h*(0.5*(obj.X2V.*(obj.sang_x4.*obj.cang_x5.*obj.cang_x6 - obj.cang_x4.*obj.sang_x5.*obj.sang_x6) ...
                -obj.X1V.*(obj.cang_x4.*obj.sang_x5.*obj.cang_x6 + obj.sang_x4.*obj.cang_x5.*obj.sang_x6) ...
                +obj.X3V.*x7));
            
            %x4_next
            x7 = x7 + obj.h*(0.5*(-obj.X1V.*(obj.sang_x4.*obj.cang_x5.*obj.cang_x6 - obj.cang_x4.*obj.sang_x5.*obj.sang_x6) ...
                -obj.X2V.*(obj.cang_x4.*obj.sang_x5.*obj.cang_x6 + obj.sang_x4.*obj.cang_x5.*obj.sang_x6) ...
                -obj.X3V.*(obj.cang_x4.*obj.cang_x5.*obj.sang_x6 - obj.sang_x4.*obj.sang_x5.*obj.cang_x6) ));
            %
            %
            size_X5 = size(x7); %will be used in reshaping yaw, pitch, roll new to original n-D matrix
            
            obj.X4_next = obj.X4_next(:);
            obj.X5_next = obj.X5_next(:);
            obj.X6_next = obj.X6_next(:);
            x7 = x7(:);
            
            Qsquared_sum = sqrt( obj.X4_next.^2 + obj.X5_next.^2 + obj.X6_next.^2 + x7.^2 );
            
            
            obj.X4_next = obj.X4_next./Qsquared_sum;
            obj.X5_next = obj.X5_next./Qsquared_sum;
            obj.X6_next = obj.X6_next./Qsquared_sum;
            x7 = x7./Qsquared_sum;
            
            r1_yaw_next = atan2( 2.*(obj.X6_next.*obj.X5_next + x7.*obj.X4_next),...
                x7.^2 + obj.X6_next.^2 - obj.X5_next.^2 - obj.X4_next.^2 );
            r2_pitch_next = asin( -2.*(obj.X6_next.*obj.X4_next - x7.*obj.X5_next) );
            r3_roll_next = atan2(2.*(obj.X5_next.*obj.X4_next + x7.*obj.X6_next), ...
                x7.^2 - obj.X6_next.^2 - obj.X5_next.^2 + obj.X4_next.^2);
            
            obj.X4_next = reshape(r1_yaw_next, size_X5);
            obj.X5_next = reshape(r2_pitch_next, size_X5);
            obj.X6_next = reshape(r3_roll_next, size_X5);
            %temporary repmat%%%%%%%%%%%%%%%%%%%%%%%%%%%
            nq = obj.n_mesh_q;
            %             nw = obj.n_mesh_w;
            nu = length(obj.U_vector);
            
            obj.X1_next = repmat(obj.X1_next,[1 1 1 nq nq nq 1 nu nu]); %X1_next size = nw x nw x nw x 1 x 1 x 1 x nu
            obj.X2_next = repmat(obj.X2_next,[1 1 1 nq nq nq nu 1 nu]); %X2_next size = nw x nw x nw x 1 x 1 x 1 x 1 x nu
            obj.X3_next = repmat(obj.X3_next,[1 1 1 nq nq nq nu nu 1]); %X3_next size = nw x nw x nw x 1 x 1 x 1 x 1 x 1 x nu
            obj.X4_next = repmat(obj.X4_next,[1 1 1 1 1 1 nu nu nu]); %X4_next size = nw x nw x nw x nq x nq x nq
            obj.X5_next = repmat(obj.X5_next,[1 1 1 1 1 1 nu nu nu]); %X5_next size = nw x nw x nw x nq x nq x nq
            obj.X6_next = repmat(obj.X6_next,[1 1 1 1 1 1 nu nu nu]); %X6_next size = nw x nw x nw x nq x nq x nq
            %%%%%%%%%%%%
        end
        
        function linear_control_response(spacecraft, X0, T_final, dt)
            %example:  S.linear_control_response(S.defaultX0, 100,1e-3
            
            if nargin < 2
                %sample initial state
                X0 = spacecraft.defaultX0;
                T_final = spacecraft.T_final ;
                dt = spacecraft.h;
            end
            
            N = T_final/dt;
            U = spacecraft.U_vector';
            X(:,1) = X0;
            qc = [1, 0, 0, 0;...
                0, 1, 0, 0;...
                0, 0, 1, 0;...
                0, 0, 0, 1]; % q command (at origin, is equal to I(4x4) )
            K = [0.2, 0, 0;...
                0, 0.2, 0;
                0, 0, 0.2];
            C = [1, 0, 0;...
                0, 1, 0;
                0, 0, 1];
            %             keyboard;
            tic
            for k_stage=1:N
                % qe = qc*q;
                q = X(4:7, k_stage);
                qe = qc*q;
                w = X(1:3, k_stage);
                U(:,k_stage) = -K*qe(1:3) - C*w;
                X(:,k_stage+1) = next_stage_states(spacecraft, [w',q'], U(:,k_stage)', dt);
                [x_yaw,x_pitch,x_roll] = quat2angle([X(7,k_stage),X(6,k_stage),X(5,k_stage),X(4,k_stage)]);
                XANGLES(:,k_stage) = [x_yaw;x_pitch;x_roll];
            end
            q_squared_sum = sqrt(X(4,:).^2 + X(5,:).^2 + X(6,:).^2 + X(7,:).^2); %check quaternions
            %print time and error
            % note: quaternions deviation from (sum(Q.^2) = 1) at T_final is a measure of error in ode solver
            fprintf(...
                'Done - Time elapsed for caculations: %f - states max error: %.5g\n',...
                toc, sqrt((q_squared_sum(end) - 1)))
            
            time_v = linspace(0, T_final, N); %plot time vector
            % plot controls
            figure
            plot(time_v, U(1,:),'--')
            hold on
            plot(time_v, U(2,:),'--')
            plot(time_v, U(3,:),'--')
            grid on
            title('Control Inputs')
            legend('u1','u2','u3')
            
            %plot states
            figure
            hold on
            grid on
            for n_state = [1 2 3]
                plot(time_v, X(n_state, 1:end-1).*180/pi)
            end
            legend('\omega_1','\omega_2','\omega_3')
            xlabel('time (s)')
            ylabel('\omega (deg/sec)')
            
            figure
            hold on
            grid on
            for n_state = [4 5 6 7]
                plot(time_v, X(n_state, 1:end-1))
            end
            legend('q1','q2','q3','q4')
            xlabel('time (s)')
            
            figure
            hold on
            grid on
            for n_state = [1 2 3]
                plot(time_v, XANGLES(n_state, :).*180/pi)
            end
            legend('\theta_1','\theta_2','\theta_3')
            xlabel('time (s)')
            ylabel('\theta (deg)')
            
        end
        
        function [X1_dot,X2_dot,X3_dot,...
                X5_dot,X6_dot,X7_dot] = spacecraft_dynamics(obj,x1,x2,x3,x4,x5,x6,u1,u2,u3)
            %returns the derivatives x_dot = f(X,u)
            %J is assumed to be diagonal, J12 = J23 = ... = 0
            
        end
        
        function X_dot = spacecraft_dynamics_list(obj, X,U)
            %returns the derivatives x_dot = f(X,u)
            % FIXED q~ convention
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
            X_dot(:,4) = 0.5*(x3.*x5 -x2.*x6 +x1.*x7);
            X_dot(:,5) = 0.5*(-x3.*x4 +x1.*x6 +x2.*x7);
            X_dot(:,6) = 0.5*(x2.*x4 -x1.*x5 +x3.*x7);
            X_dot(:,7) = 0.5*(-x1.*x4 -x2.*x5 -x3.*x6);
        end
        
        %% functions for simplified dynamics
        function [w_next,t_next] = next_stage_states_simplified(obj,W, T, U, J, h)
            w_next = RK4_w(obj,W, U, J, h);
            t_next = RK4_t(obj,T, W, h);
        end
        
        function W2 = RK4_w(obj,w, U, J, h)
            %calculates next stage (k+1) states
            % X2 = X1 + dt*a_d(X,u) where a_d is the spacecraft dynamics
            %first order taylor expansion
            %X2 = X1 + dt*spacecraft_dynamics(spacecraft, X1, U);
            
            % Runge-Kutta - 4th order
            % h = dt;
            k1 = wdynamics(obj,w , U, J);
            k2 = wdynamics(obj,(w + k1*h/2), U, J);
            k3 = wdynamics(obj,(w + k2*h/2), U, J);
            k4 = wdynamics(obj,(w + k3*h), U, J);
            
            W2 = w + h*(k1 + 2*k2 + 2*k3 + k4)/6;
        end
        
        function T2 = RK4_t(obj,T1, W1, h)
            %calculates next stage (k+1) states
            % X2 = X1 + dt*a_d(X,u) where a_d is the spacecraft dynamics
            %first order taylor expansion
            %X2 = X1 + dt*spacecraft_dynamics(spacecraft, X1, U);
            
            % Runge-Kutta - 4th order
            % h = dt;
            k1 = tdynamics(obj,W1);
            k2 = tdynamics(obj,(W1 + k1*h/2));
            k3 = tdynamics(obj,(W1 + k2*h/2));
            k4 = tdynamics(obj,(W1 + k3*h));
            
            T2 = T1 + h*(k1 + 2*k2 + 2*k3 + k4)/6;
        end
        
        function w_dot = wdynamics(~,~, U, J)
            w_dot = U/J;
        end
        function t_dot = tdynamics(~,w)
            t_dot = w;
        end
        %% end of simplified dynamics functions
        
        function X2 = next_stage_states(spacecraft, X1, U, h, mode)
            if nargin < 5
                mode = 'RK4';
            end
            %calculates next stage (k+1) states
            % X2 = X1 + dt*a_d(X,u) where a_d is the spacecraft dynamics
            %first order taylor expansion
            %X2 = X1 + dt*spacecraft_dynamics(spacecraft, X1, U);
            
            % Runge-Kutta - 4th order
            % h = dt;
            k1 = spacecraft_dynamics_list(spacecraft, X1 , U);
            if (strcmp(mode,'RK4'))
                k2 = spacecraft_dynamics_list(spacecraft, (X1 + k1*h/2), U);
                k3 = spacecraft_dynamics_list(spacecraft, (X1 + k2*h/2), U);
                k4 = spacecraft_dynamics_list(spacecraft, (X1 + k3*h), U);
                
                
                X2 = X1 + h*(k1 + 2*k2 + 2*k3 + k4)/6;
                
            elseif(strcmp(mode,'taylor'))
                X2 = X1 + h*k1;
            end
            q_sqsum = sqrt(X2(4).^2 +  X2(5).^2 +  X2(6).^2 +  X2(7).^2);
            X2(4:7) = X2(4:7)./q_sqsum;
            
        end
        
        function [sr_5, sr_6, sr_7] = mesh_quaternion(this)
            obj.s_yaw = linspace(deg2rad(this.yaw_min), deg2rad(this.yaw_max), this.n_mesh_q);
            obj.s_pitch = linspace(deg2rad(this.pitch_min), deg2rad(this.pitch_max), this.n_mesh_q);
            obj.s_roll = linspace(deg2rad(this.roll_min), deg2rad(this.roll_max), this.n_mesh_q);
            
            angles = [s_yaw(:) s_pitch(:) s_roll(:)];
            
            cang = cos( angles/2 );
            sang = sin( angles/2 );
            
            %             sr_5 = cang(:,1).*sang(:,2).*cang(:,3) + sang(:,1).*cang(:,2).*sang(:,3);
            %             sr_6 = cang(:,1).*cang(:,2).*sang(:,3) - sang(:,1).*sang(:,2).*cang(:,3);
            %             sr_7 = cang(:,1).*cang(:,2).*cang(:,3) + sang(:,1).*sang(:,2).*sang(:,3);
            sr_5 = sr_5';
            sr_6 = sr_6';
            sr_7 = sr_7';
            
        end
        
        function reshape_states(obj)
            n_mesh_w1 = length(obj.sr_1);
            n_mesh_w2 = length(obj.sr_2);
            n_mesh_w3 = length(obj.sr_3);
            
            obj.X1V = reshape(obj.sr_1, [n_mesh_w1 1]);
            obj.X2V = reshape(obj.sr_2, [1 n_mesh_w2]);
            obj.X3V = reshape(obj.sr_3, [1 1 n_mesh_w3]);
            
            n_mesh_yaw = length(obj.s_yaw);
            obj.cang_x4 = reshape( cos( obj.s_yaw/2 ), [1 1 1 n_mesh_yaw 1 1]);
            obj.sang_x4 = reshape( sin( obj.s_yaw/2 ), [1 1 1 n_mesh_yaw 1 1]);
            
            n_mesh_pitch = length(obj.s_pitch);
            obj.cang_x5 = reshape( cos( obj.s_pitch/2 ), [1 1 1 1 n_mesh_pitch 1]);
            obj.sang_x5 = reshape( sin( obj.s_pitch/2 ), [1 1 1 1 n_mesh_pitch 1]);
            
            n_mesh_roll = length(obj.s_roll);
            obj.cang_x6 = reshape( cos( obj.s_roll/2 ), [1 1 1 1 1 n_mesh_roll]);
            obj.sang_x6 = reshape( sin( obj.s_roll/2 ), [1 1 1 1 1 n_mesh_roll]);
            
            n_mesh_u = length(obj.U_vector);
            obj.U1V = reshape(obj.U_vector, [1 1 1 1 1 1 n_mesh_u]);
            obj.U2V = reshape(obj.U_vector, [1 1 1 1 1 1 1 n_mesh_u]);
            obj.U3V = reshape(obj.U_vector, [1 1 1 1 1 1 1 1 n_mesh_u]);
        end
        
        function get_optimal_path(obj, X0, method)
            if nargin < 2
                X0 = obj.defaultX0;
                method = 'nearest';
            end
            X = zeros(7,obj.N_stage);
            X(:,1) = X0;
            U = zeros(3,obj.N_stage);
            X_ANGLES = zeros(9,obj.N_stage);
            tic
            for k_stage=1:obj.N_stage-1
                % qe = qc*q;
                %% proper q input to quat2angle()
                [x_yaw,x_pitch,x_roll] = quat2angle([X(7,k_stage),X(6,k_stage),X(5,k_stage),X(4,k_stage)]);
                %                 X(4:7,k_stage) = angle2quat(x_yaw,x_pitch,x_roll)';
                FU1 = griddedInterpolant({obj.sr_1, obj.sr_2, obj.sr_3, obj.s_yaw, obj.s_pitch, obj.s_roll},...
                    obj.U1_Opt,method);
                FU2 = griddedInterpolant({obj.sr_1, obj.sr_2, obj.sr_3, obj.s_yaw, obj.s_pitch, obj.s_roll},...
                    obj.U2_Opt,method);
                FU3 = griddedInterpolant({obj.sr_1, obj.sr_2, obj.sr_3, obj.s_yaw, obj.s_pitch, obj.s_roll},...
                    obj.U3_Opt,method);
                
                U(1,k_stage) = FU1(...
                    X(1,k_stage),X(2,k_stage),X(3,k_stage), ...
                    x_yaw,x_pitch,x_roll);
                
                U(2,k_stage) = FU2(...
                    X(1,k_stage),X(2,k_stage),X(3,k_stage), ...
                    x_yaw,x_pitch,x_roll);
                
                U(3,k_stage) = FU3(...
                    X(1,k_stage),X(2,k_stage),X(3,k_stage), ...
                    x_yaw,x_pitch,x_roll);
                
                %uses taylor for next stage calculations
                X(:,k_stage+1) = next_stage_states(obj, X(:,k_stage)', U(:,k_stage)', obj.h, 'taylor');
                
                X_ANGLES(:,k_stage) = [X(1,k_stage);X(2,k_stage);X(3,k_stage);...
                    rad2deg(x_roll);rad2deg(x_pitch);rad2deg(x_yaw);...
                    U(:,k_stage)];
            end
            toc
            %% check vector values (should start from 0)
            
            v_plot = 0:obj.h:obj.T_final-obj.h;
            figure
            hold on
            for i=1:3
                plot(v_plot, U(i,:),'--')
            end
            legend('u1','u2','u3')
            grid on
            xlabel('time (s)')
            
            %plot w1,w2,w3
            figure
            hold on
            grid on
            for n_state = [1 2 3]
                plot(v_plot, X(n_state, :).*180/pi)
            end
            legend('\omega_1','\omega_2','\omega_3')
            xlabel('time (s)')
            ylabel('\omega (deg/sec)')
            
            
            %plot angles
            figure
            hold on
            for i=4:6
                plot(v_plot, X_ANGLES(i,:))
            end
            title('rotation angles')
            legend('roll','pitch','yaw')
            grid on
            ylabel('deg')
            xlabel('time (s)')
            
            
            figure
            hold on
            for i=1:6
                plot(v_plot, X(i,1:end))
            end
            title('states')
            legend('x1','x2','x3','x5','x6','x7')
            grid on
            xlabel('time (s)')
            
        end
        
        function get_optimal_path_simplified_testode45(obj)
            tic
            X0 = obj.defaultX0;
            FU1 = obj.U1_Opt;
            FU2 = obj.U2_Opt;
            FU3 = obj.U3_Opt;
            tspan = 0:obj.h:obj.T_final;
            X_ode45 = zeros(obj.N_stage, 7);
            X_ode45(1,:) = X0;
            for k_stage=1:obj.N_stage-1
                %determine U
                X_stage = X_ode45(k_stage,:);
                U1(1) = FU1(X_stage(1), 2*asin(X_stage(4)));
                U1(2) = FU2(X_stage(2), 2*asin(X_stage(5)));
                U1(3) = FU3(X_stage(3), 2*asin(X_stage(6)));
                %
                [~,X_temp] = ode45(@ode_eq,[tspan(k_stage), tspan(k_stage+1)], X_ode45(k_stage,:));
                X_ode45(k_stage+1,:) = X_temp(end,:);
            end
            toc
            T_ode45 = tspan(1:end-1)';
            [t_yaw,t_pitch,t_roll] = quat2angle([X_ode45(:,7),X_ode45(:,6),X_ode45(:,5),X_ode45(:,4)]);
            theta1_ode45 = t_yaw;
            theta2_ode45 = t_pitch;
            theta3_ode45 = t_roll;
            U1_ode45 = FU1(X_ode45(:,1), theta1_ode45);
            U2_ode45 = FU2(X_ode45(:,2), theta2_ode45);
            U3_ode45 = FU3(X_ode45(:,3), theta3_ode45);
            
            %% plot rotational speeds
            figure
            hold on
            grid on
            for n_state = [1 2 3]
                plot(T_ode45, X_ode45(:, n_state).*180/pi)
            end
            legend('\omega_1','\omega_2','\omega_3')
            xlabel('time (s)')
            ylabel('\omega (deg/sec)')
            
            %% plot angles
            figure
            hold on
            plot(T_ode45,theta1_ode45)
            plot(T_ode45,theta2_ode45)
            plot(T_ode45,theta3_ode45)
            
            title('rotation angles')
            legend('\theta_1','\theta_2','\theta_3')
            grid on
            xlabel('time (s)')
            
            %% plot controls
            figure
            plot(T_ode45, U1_ode45,'--')
            hold on
            plot(T_ode45, U2_ode45,'--')
            plot(T_ode45, U3_ode45,'--')
            grid on
            title('Control Inputs')
            legend('u1','u2','u3')
            
            %function declarations
            function x_dot = ode_eq(~,X1)
                
                x_dot = system_dynamics(X1, U1');
                x_dot = x_dot';
                function X_dot = system_dynamics(X, U)
                    x1 = X(1);
                    x2 = X(2);
                    x3 = X(3);
                    x4 = X(4);
                    x5 = X(5);
                    x6 = X(6);
                    x7 = X(7);
                    %InertiaM is a complete matrix
                    w = [x1;x2;x3];
                    w_dot = obj.InertiaM\(U - cross(w, obj.InertiaM*w));
                    
                    X_dot(1) = w_dot(1);
                    X_dot(2) = w_dot(2);
                    X_dot(3) = w_dot(3);
                    
                    X_dot(4) = 0.5*(x3.*x5 -x2.*x6 +x1.*x7);
                    X_dot(5) = 0.5*(-x3.*x4 +x1.*x6 +x2.*x7);
                    X_dot(6) = 0.5*(x2.*x4 -x1.*x5 +x3.*x7);
                    X_dot(7) = 0.5*(-x1.*x4 -x2.*x5 -x3.*x6);
                end
                
            end
        end
        
        function clear_up(obj)
            % clears up RAM from big matrices that were used in calculations and are no longer needed
            obj.X2_next = [];
            obj.X1_next = [];
            obj.X3_next = [];
            obj.X4_next = [];
            obj.X5_next = [];
            obj.X6_next = [];
            obj.J_current_state_fix = [];
            obj.F = [];
        end
    end
    
end

