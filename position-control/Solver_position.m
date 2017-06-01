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
                
                this.v_min = -10;
                this.v_max = +10;
              
                this.x_min = -15;
                this.x_max = 15;
                
                this.n_mesh_v = 100;
                this.n_mesh_x = 300;
                
                this.Mass = 4.16;
                
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
            this.U_vector = [-0.13 0 0.13]; % Thruster Force
            
            %u_opt
            size_Umat = [this.n_mesh_x, this.n_mesh_v];
            this.U1_Opt = zeros(size_Umat,'single');
            this.U2_Opt = zeros(size_Umat,'single');
            this.U3_Opt = zeros(size_Umat,'single');
            
        end
        
        function simplified_run(obj)
            
            %% mesh generation
            s_x1 = linspace(obj.x_min, obj.x_max, obj.n_mesh_x);
            s_x2 = linspace(obj.x_min, obj.x_max, obj.n_mesh_x);
            s_x3 = linspace(obj.x_min, obj.x_max, obj.n_mesh_x);
            
            s_v1 = linspace(obj.v_min, obj.v_max, obj.n_mesh_v);
            s_v2 = linspace(obj.v_min, obj.v_max, obj.n_mesh_v);
            s_v3 = linspace(obj.v_min, obj.v_max, obj.n_mesh_v);
            
            
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
        
    end
end


