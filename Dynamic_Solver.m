classdef Dynamic_Solver
    %DYNAMIC_SOLVER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        A
        B
        J
        H
        R
        Q
        N %Number of stages
        S %Number of state values
        C %Number of control values
        X
        U
        x_min
        x_max
        u_min
        u_max
        dx
        du
        u_star
        J_star
    end
    
    methods
        
        function obj = Dynamic_Solver()
            obj.Q = [0.25, 0; 0, 0.25];
            obj.A = [0.9974, 0.0539; -0.1078, 1.1591];
            obj.B = [0.0013; 0.0539];
            obj.R = 0.05;
            obj.N = 5;
            obj.S = 2;
            obj.C = 1;
            obj.X = zeros(obj.S,obj.N);
            obj.U = zeros(obj.C,obj.N);
            obj.dx = 100;
            obj.du = 100;
            obj.x_max = 8;
            obj.x_min = -8;
            obj.u_max = 8;
            obj.u_min = -12;
            obj.X(:,obj.N) = [0; 0];

            
        end
        
        function obj = run(obj)
            % K = 0
            % Calculate and store J*NN = h(xi(N)) for all x(N)
            s_r = linspace(obj.x_min,obj.x_max,obj.dx);
            c_r = linspace(obj.u_min,obj.u_max,obj.du);
            
            [X1_mesh, X2_mesh] = ndgrid(s_r, s_r);
            U_mesh = linspace(obj.u_max, obj.u_min,obj.du);
            obj.J_star = zeros([size(X1_mesh),obj.N]);
            obj.u_star = obj.J_star;
            % Increase K by 1
            for k=1:obj.N-1
                tic
                for i=1:obj.dx % Set xi(N-k) == starting quantized value by making i = 1
                    % set COSMIN to a large positive number
                    X1 = X1_mesh(i);
                    
                    for ii=obj.dx
                    X2 = X2_mesh(ii);
                    
                    COSTMIN = 1000000; %set to a finite large number to increase performance
                    UMIN = [];
                    % set ui(N-k) to the starting quantized value by making j = 1
                    for j=1:obj.du
                        Ui = U_mesh(j);
                        % Calculate the value of x(i,j)(N-k +1) = a_D(xi(N -k),uj(N-k))
                        X_next = obj.A*[X1;X2] + obj.B*Ui;
                        % Use this value of x(i,j)(N-k+1) to select the appropriate
                        % stored value of J*{(N-k),N} (x(i,j)(N-k+1))
                        % if x(i,j)(N-k+1) is not a grid value, interpolation is
                        % required.
                        %
                        % using function $interpn$
                        %J_opt_next = interpn(X1_mesh, X2_mesh,obj.J_star(:,:,obj.N-k+1),...
                        %    X_next(1),X_next(2));
                        %
                        %using gridded
                        F = griddedInterpolant(X1_mesh, X2_mesh,...
                                                obj.J_star(:,:,obj.N-k+1),'linear');
                        J_opt_next = F(X_next(1),X_next(2));
                        
                        %
                        % Compute C*(N-k,N) (xi(N-k),uj(N-k)) the minimum cost over the
                        % final (N-k) stages if uj(N-k) is applied at xj(N-k)
                        %
                        %C_star = obj.g_D(X1_mesh(i),X2_mesh(i),U_mesh(j),obj.Q,obj.R) ... +
                        %         + J_opt_next ;
                        %
                         C_star = [X1;X2]' * obj.Q * [X1;X2] ...
                             + Ui' * obj.R * Ui + J_opt_next ;
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
                    obj.u_star(i,ii,obj.N-k) = UMIN;
                    % store COSMIN in COST(N-k,i)
                    obj.J_star(i,ii,obj.N-k) = COSTMIN;
                    
                    end  %end of for loop for X2 1:dx
                end    %end of for loop when i = S
                fprintf('step %d - %f seconds\n', k, toc)
            end %end of for loop when k = N
                      
        end
        
        function X1_new = a_D(X1,X2,Ui)
            keyboard;
             X1_new = obj.A*[X1;X2] + obj.B*Ui;
%             X1_new = [A(1).*X1+ A(3).*X2 + B(1).*Ui; A(2).*X1+ A(4).*X2 + B(2).*Ui];
        end
        
        function J = g_D(X1,X2,Ui,Q,R)
           J = [X1;X2]' * Q * [X1;X2] + Ui' * R * Ui;
        end
        
        function get_optimal_path(X)
            %To be Completed
            %Print Optimal controls, UOPT(N-k,I) and min costs, COST(N-k,I)
            %for all quantized state points (I = 1,2,..,S) and all stages
            % K = 1:N
        end
    end
end

