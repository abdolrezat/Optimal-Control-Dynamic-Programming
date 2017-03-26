classdef Dynamic_Solver
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
            obj.dx = 35;
            obj.du = 100;
            obj.x_max = 3;
            obj.x_min = -2.5;
            obj.u_max = 10;
            obj.u_min = -40;
        end
        
        function obj = run(obj)
            try
            % K = 0
            % Calculate and store J*NN = h(xi(N)) for all x(N)
            s_r = linspace(obj.x_min,obj.x_max,obj.dx);
            [obj.X1_mesh, obj.X2_mesh] = ndgrid(s_r, s_r);
            
            U_mesh = linspace(obj.u_min, obj.u_max, obj.du);
            obj.J_star = zeros([size(obj.X1_mesh),obj.N]);
            obj.u_star = obj.J_star;
                        % Increase K by 1
            for k=1:obj.N-1
                tic
                F = griddedInterpolant(obj.X1_mesh, obj.X2_mesh,...
                                                obj.J_star(:,:,obj.N-k+1),'linear');
                for i1=1:obj.dx % Set xi(N-k) == starting quantized value by making i = 1
                    for i2=1:obj.dx
                    X1 = obj.X1_mesh(i1,i2);
                    X2 = obj.X2_mesh(i1,i2);
                    
                    % set COSMIN to a large positive number
                    COSTMIN = 10000; %set to a finite large number to increase performance
                    UMIN = [];
                    % set ui(N-k) to the starting quantized value by making j = 1
                    
                    for jj=1:obj.du
                        Ui = U_mesh(jj);
                        
                        % Calculate the value of x(i,j)(N-k +1) = a_D(xi(N -k),uj(N-k))
                        X_next = a_D(obj,X1,X2,Ui);
                        % Use this value of x(i,j)(N-k+1) to select the appropriate
                        % stored value of J*{(N-k),N} (x(i,j)(N-k+1))
                        % if x(i,j)(N-k+1) is not a grid value, interpolation is
                        % required.
                        %
                        % using function $interpn$
                        %J_opt_next = interpn(X1_mesh, X2_mesh,obj.J_star(:,:,obj.N-k+1),...
                        %    X_next(1),X_next(2));
                        %
                        %using griddedInterpolant
                        
                        J_opt_next = F(X_next(1),X_next(2));
                        
                        %
                        % Compute C*(N-k,N) (xi(N-k),uj(N-k)) the minimum cost over the
                        % final (N-k) stages if uj(N-k) is applied at xj(N-k)
                        %
                        %C_star = obj.g_D(X1_mesh(i),X2_mesh(i),U_mesh(j),obj.Q,obj.R) ... +
                        %         + J_opt_next ;
                        %
                        % C_star = [X1;X2]' * obj.Q * [X1;X2] ...
                        %     + Ui' * obj.R * Ui + J_opt_next ;
                        
                        C_star = g_D(obj,X1,X2,Ui) + J_opt_next ;

                        % if C* just calculated less than COSMIN store this value as
                        % COSMIN and store the value uj(N-k) in UMIN
                        if(C_star < COSTMIN)
                            COSTMIN = C_star;
                            UMIN = Ui;
                            
                            % debug start 
                            %-- check u limit -- disable for performance
                            %improvement
                            %if(UMIN == obj.u_max) 
                            %    warning('U control reached limit x1 = %f, x2 = %f, k = %d', X1, X2, k);
                            %end
                            %-- save X_next records
                            %X_N(i,ii) = sum(X_next);
                            % debug end
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
            
            catch e
                e.throw
            end        
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
                keyboard
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
        
                   
        function X1_new = a_D(obj,X1,X2,Ui)
            %keyboard;
             X1_new = obj.A*[X1;X2] + obj.B*Ui;
%             X1_new = [A(1).*X1+ A(3).*X2 + B(1).*Ui; A(2).*X1+ A(4).*X2 + B(2).*Ui];
        end
        
        function J = g_D(obj,X1,X2,Ui)
           J = [X1;X2]' * obj.Q * [X1;X2] + Ui' * obj.R * Ui;
        end
        
    end
    
    methods (Static)

        function b = compare_data(obj1,obj2)
           % use this function to compare saved datas 
           % check J* matrix
           if( isempty(obj1.J_star) || isempty(obj2.J_star) )
               error('stop throwing empty data at me')
           end
           
           c1 = obj1.J_star == obj2.J_star ;
           if(all(c1(:)))
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

