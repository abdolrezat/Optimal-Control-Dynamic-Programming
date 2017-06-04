classdef Solver_position < handle
    %SOLVER_POSITION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        N % number of stages
        Mass % Mass
        
        size_state_mat
        
        U_vector % values that are applied as control inputs
        % big data
        %J_next_states_opt single  % optimum values of J will be stored in this matrix and recalled each stage
        J_current_state_fix single  % values of cost to reach next stage, always fixed in Time Invariant system
        
        v_min
        v_max
        
        x_min
        x_max
        n_mesh_v
        n_mesh_x
        Qx1
        Qx2
        Qx3
        Qv1
        Qv2
        Qv3
        
        R1
        R2
        R3
        
        T_final
        h
        N_stage
        defaultX0
        
        U1_Opt
        U2_Opt
        U3_Opt
        
    end
    
    methods
        function this = Solver_position()
            if nargin < 1
                
                this.v_min = -0.5;
                this.v_max = +0.5;
                
                this.x_min = -.5;
                this.x_max = .5;
                
                this.n_mesh_v = 200;
                this.n_mesh_x = 200;
                
                this.Mass = 4.16;
                
                this.Qx1 = 6;
                this.Qx2 = 6;
                this.Qx3 = 6;
                this.Qv1 = 6;
                this.Qv2 = 6;
                this.Qv3 = 6;
                
                this.R1 = 0.1;
                this.R2 = 0.1;
                this.R3 = 0.1;
                
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
            this.U_vector = [-0.13 0 0.13]*2; % Thruster Force
            
            %u_opt
            size_Umat = [this.n_mesh_x, this.n_mesh_v];
            this.U1_Opt = zeros(size_Umat,'single');
            this.U2_Opt = zeros(size_Umat,'single');
            this.U3_Opt = zeros(size_Umat,'single');
            
        end
        
        function simplified_run(obj)
            
            %% mesh generation
            s_x1 = sym_linspace(obj, obj.x_min, obj.x_max, obj.n_mesh_x);
            s_x2 = sym_linspace(obj, obj.x_min, obj.x_max, obj.n_mesh_x);
            s_x3 = sym_linspace(obj, obj.x_min, obj.x_max, obj.n_mesh_x);
            obj.n_mesh_x = length(s_x1);
            s_v1 = sym_linspace(obj, obj.v_min, obj.v_max, obj.n_mesh_v);
            s_v2 = sym_linspace(obj, obj.v_min, obj.v_max, obj.n_mesh_v);
            s_v3 = sym_linspace(obj, obj.v_min, obj.v_max, obj.n_mesh_v);
            obj.n_mesh_v = length(s_v1);
            
            
            %% initialization
            %states, controls
            [X1,V1,U1] = ndgrid(s_x1,s_v1,obj.U_vector);
            [X2,V2,U2] = ndgrid(s_x2,s_v2,obj.U_vector);
            [X3,V3,U3] = ndgrid(s_x3,s_v3,obj.U_vector);
            
            J_current = @(x,v,U,Qx,Qv,R)(Qx * x.^2 + Qv * v.^2 + R*U.^2);
            %calculating J fixed
            J_current_1 = J_current(X1,V1,U1,obj.Qx1,obj.Qv1,obj.R1);
            F1 = griddedInterpolant({s_x1,s_v1},zeros(obj.n_mesh_x,obj.n_mesh_v),'linear');
            
            J_current_2 = J_current(X2,V2,U2,obj.Qx2,obj.Qv2,obj.R2);
            F2 = griddedInterpolant({s_x2,s_v2},zeros(obj.n_mesh_x,obj.n_mesh_v),'linear');
            
            J_current_3 = J_current(X3,V3,U3,obj.Qx3,obj.Qv3,obj.R3);
            F3 = griddedInterpolant({s_x3,s_v3},zeros(obj.n_mesh_x,obj.n_mesh_v),'linear');
            
            %calculating next stage states
            %
            [x1_next,v1_next] = next_stage_states_simplified(obj, X1, V1, U1, obj.h);
            [x2_next,v2_next] = next_stage_states_simplified(obj, X2, V2, U2, obj.h);
            [x3_next,v3_next] = next_stage_states_simplified(obj, X3, V3, U3, obj.h);
            
            %beginning (reverse) stage calculations
            whandle = waitbar(0,'Calculation in Progress...');
            for k_s = obj.N_stage-1:-1:1
                tic
                %     C = F1(w1_next,t1_next);
                [F1.Values, U1_idx] = min( J_current_1 + F1(x1_next,v1_next), [], 3);
                [F2.Values, U2_idx] = min( J_current_2 + F2(x2_next,v2_next), [], 3);
                [F3.Values, U3_idx] = min( J_current_3 + F3(x3_next,v3_next), [], 3);
                
                waitbar( 1 - k_s/obj.N_stage, whandle);
                fprintf('step %d - %f seconds\n', k_s, toc)
            end
            
            %set U* Optimal, replace Fi itself to save memory
            obj.U1_Opt = griddedInterpolant({s_x1,s_v1}, obj.U_vector(U1_idx),'nearest');
            obj.U2_Opt = griddedInterpolant({s_x2,s_v2}, obj.U_vector(U2_idx),'nearest');
            obj.U3_Opt = griddedInterpolant({s_x3,s_v3}, obj.U_vector(U3_idx),'nearest');
            close(whandle)
            fprintf('stage calculation complete!\n')
            
        end
        
        function [x_next,v_next] = next_stage_states_simplified(obj, X, V, U, h)
            x_next = RK4_x(obj, X, V, h);
            v_next = RK4_v(obj, V, U, h);
        end
        
        function X2 = RK4_x(obj, X1, V, h)
            % Runge-Kutta - 4th order
            % h = dt;
            k1 = xdynamics(obj, V);
            k2 = xdynamics(obj,(V + k1*h/2));
            k3 = xdynamics(obj,(V + k2*h/2));
            k4 = xdynamics(obj,(V + k3*h));
            
            X2 = X1 + h*(k1 + 2*k2 + 2*k3 + k4)/6;
            
        end
        
        function x_dot = xdynamics(~,v)
            x_dot = v;
        end
        
        function V2 = RK4_v(obj, V1, U, h) % does not need RK4, ki's are equal
            % Runge-Kutta - 4th order
            % h = dt;
            k1 = vdynamics(obj, V1 , U);
            k2 = vdynamics(obj,(V1 + k1*h/2), U);
            k3 = vdynamics(obj,(V1 + k2*h/2), U);
            k4 = vdynamics(obj,(V1 + k3*h), U);
            
            V2 = V1 + h*(k1 + 2*k2 + 2*k3 + k4)/6;
        end
        
        function v_dot = vdynamics(obj,~, U)
            v_dot = U/obj.Mass;
        end
        
        
        function get_optimal_path(obj)
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
                [R,V] = update_RV_target(obj, R0, V0, tspan(k_stage));
                norm_R = (R*R')^.5; %norm R
                RdotV = sum(R.*V); %dot product
                crossRV = [R(2).*V(3)-R(3).*V(2); % cross product of R and V
                    R(3).*V(1)-R(1).*V(3);
                    R(1).*V(2)-R(2).*V(1)];
                H  = (crossRV'*crossRV)^.5 ; %norm(crossRV);
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
        
        function v = sym_linspace(~,a,b,n)
            if(a>0)
                error('minimum states are not negative, use normal linspace')
            end
            v_1 = linspace(a,0,ceil((n)/2)+1);
            v_2 = linspace(0,b,ceil((n)/2)+1);
            v_2 = v_2(2:end); %remove first zero
            v = [v_1,v_2];
        end
    end
    
    methods(Static)
        
    end
    
end


