classdef Dynamic_Solver < handle
    %DYNAMIC_SOLVER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        A
        B
        H
        R
        Q
        N %Number of stages
        S %Number of state values
        C %Number of control values
        %X
        %U
        x_min
        x_max
        u_min
        u_max
        dx
        du
        u_star
        J_star
        X1_mesh
        X2_mesh
    end
    
    methods
        
        function obj = Dynamic_Solver()
            obj.Q = [0.25, 0; 0, 0.05];
            obj.A = [0.9974, 0.0539; -0.1078, 1.1591];
            obj.B = [0.0013; 0.0539];
            obj.R = 0.05;
            obj.N = 130;
            obj.S = 2;
            obj.C = 1;
            %obj.X = zeros(obj.S,obj.N);
            %obj.U = zeros(obj.C,obj.N);
            obj.dx = 100;
            obj.du = 100;
            obj.x_max = 3;
            obj.x_min = -2.5;
            obj.u_max = 10;
            obj.u_min = -40;
        end
        
        function obj = run(obj)
            % K = 0
            % Calculate and store J*NN = h(xi(N)) for all x(N)
            s_r = linspace(obj.x_min,obj.x_max,obj.dx);
            [obj.X1_mesh, obj.X2_mesh] = ndgrid(s_r, s_r);
            
            U_mesh = linspace(obj.u_min, obj.u_max, obj.du);
            
            %3d grid for voctorization in calculation of C_star_M
            [X1_mesh_3D,X2_mesh_3D,U_mesh_3D] = ndgrid(s_r,s_r,U_mesh);
            %
            obj.J_star = zeros([size(obj.X1_mesh),obj.N]);
            obj.u_star = obj.J_star;
            % Increase K by 1
            for k=1:obj.N-1
                tic
                C_star_M = obj.Q(1)*X1_mesh_3D.^2 + ...
                    obj.Q(4)*X2_mesh_3D.^2 + obj.R * U_mesh_3D.^2;
                J_opt_M = C_star_M + J_opt(obj, X1_mesh_3D, X2_mesh_3D, U_mesh_3D, k);
                for i1=1:obj.dx % Set xi(N-k) == starting quantized value by making i = 1
                    for i2=1:obj.dx
                        % set COSMIN to a large positive number
                        COSTMIN = 10000; %set to a finite large number to increase performance
                        UMIN = 10000;
                        % set ui(N-k) to the starting quantized value by making j = 1
                        for jj=1:obj.du
                            Ui = U_mesh(jj);
                            C_star = J_opt_M(i1,i2,jj) ;
                            % if C* just calculated less than COSMIN store this value as
                            % COSMIN and store the value uj(N-k) in UMIN
                            if(C_star < COSTMIN)
                                COSTMIN = C_star;
                                UMIN = Ui;
                            end
                        end%end of for loop when j = C
                        %[i_x1,i_x2] = ind2sub(size(X1_mesh),i); % used for 2D
                        %X's (e.g. X1_mesh(i){mxm}
                        % store UMIN in UOPT(N-k,I)
                        obj.u_star(i1,i2,obj.N-k) = UMIN;
                        % store COSMIN in COST(N-k,i)
                        obj.J_star(i1,i2,obj.N-k) = COSTMIN;
                        
                    end  %end of innerX2 for loop ii = 1:dx
                    
                end    %end of outerX1 for loop when i = dx
                
                fprintf('step %d - %f seconds\n', k, toc)
            end %end of for loop when k = N
            
            
        end
        
        
        function get_optimal_path(obj, X0)
            if nargin < 2
                X0 = [2; 1]
            end
            %Store Optimal controls, UOPT(N-k,I) and min costs, COST(N-k,I)
            %for all quantized state points (I = 1,2,..,S) and all stages
            % K = 1:N
            plot_k_max = obj.N;
            v = 1:plot_k_max;
            X = zeros(obj.S,obj.N);
            U = zeros(obj.C,obj.N);
            J = U;
            X(:,1) = X0;
            for k=1:obj.N-1
                
                Fu = griddedInterpolant(obj.X1_mesh, obj.X2_mesh,...
                    obj.u_star(:,:,k),'linear');
                
                U(k) = Fu(X(1,k),X(2,k));
                
                Fj = griddedInterpolant(obj.X1_mesh, obj.X2_mesh,...
                    obj.J_star(:,:,k),'linear');
                J(k) = Fj(X(1,k),X(2,k));
                
                X(:,k+1) = a_D(obj,X(1,k),X(2,k),U(k));
                % X(:,k+1) = obj.A*X(:,k) + obj.B*U(k);
            end
            k = k+1;
            %-- Optimal Control Input u*
            Fu = griddedInterpolant(obj.X1_mesh, obj.X2_mesh,...
                obj.u_star(:,:,k),'linear');
            U(k) = Fu(X(1,k),X(2,k));
            
            %-- Commented -- Cost of path
            %Fj = griddedInterpolant(obj.X1_mesh, obj.X2_mesh,...
            %   obj.J_star(:,:,k),'linear');
            %J(k) = Fj(X(1,k),X(2,k));
            
            %Print Optimal Controls
            plot(v,X(1,v))
            hold on
            plot(v,X(2,v),'r')
            plot(v,U(v),'--')
            title('Optimal control for initial state X0')
            xlabel('stage - k')
            ylabel('state and inputs')
            legend('X1', 'X2', 'u*');
            grid on
            xlim([v(1) v(end)])
            
        end
        
        
        function [Xnext_M1,Xnext_M2] = a_D_M(obj,X1,X2,Ui)
            %keyboard;
            Xnext_M1 = obj.A(1)*X1 + obj.A(3)*X2 + obj.B(1)*Ui;
            Xnext_M2 = obj.A(2)*X1 + obj.A(4)*X2 + obj.B(2)*Ui;
            %             X1_new = [A(1).*X1+ A(3).*X2 + B(1).*Ui; A(2).*X1+ A(4).*X2 + B(2).*Ui];
        end
        
        
        function X1_new = a_D(obj,X1,X2,Ui)
            X1_new = obj.A*[X1;X2] + obj.B*Ui;
        end
        
        function J = g_D(obj,X1,X2,Ui)
            J = [X1;X2]' * obj.Q * [X1;X2] + Ui' * obj.R * Ui;
        end
        
        function J_opt_M = J_opt(obj,X1,X2,Ui,k)
            F = griddedInterpolant(obj.X1_mesh, obj.X2_mesh,...
                obj.J_star(:,:,obj.N-k+1),'linear');
            %[X_next_M1,X_next_M2] = a_D_M(obj,X1_mesh_3D, X2_mesh_3D, U_mesh_3D);
            X_next_M1 = obj.A(1)*X1 + obj.A(3)*X2 + obj.B(1)*Ui;
            X_next_M2 = obj.A(2)*X1 + obj.A(4)*X2 + obj.B(2)*Ui;
            J_opt_M = F(X_next_M1,X_next_M2);
        end
    end
    
    methods (Static)
        
        function b = compare_data(obj1,obj2)
            % use this function to compare saved datas
            % check J* matrix
            if( isempty(obj1.J_star) || isempty(obj2.J_star) )
                error('stop throwing empty data at me')
            end
            %compare
            if( isequal(obj1.J_star, obj2.J_star) )
                disp('J_star matrices comparison -- Match!')
                b = true;
            else
                warning('J_star matrices -- Do NOT match')
                b = false;
            end
        end
        % end
    end
    
    
end

