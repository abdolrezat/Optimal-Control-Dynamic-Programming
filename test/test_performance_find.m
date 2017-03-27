%we'll test the performance of Two approaches: 1.using "elementwise
%multiplication + find()" to get the min of a 100x100x1000 J matrix and 2.
%using three (for) loops the usual approach
Q = [0.25, 0; 0, 0.05];
A = [0.9974, 0.0539; -0.1078, 1.1591];
B = [0.0013; 0.0539];
R = 0.05;
N = 30;

[X1, X2] = ndgrid(1:N,1:N);
U = 1:1000;
%% ----------------- Non-Vectorized Implementation ---------------------- %%

tic
for i=1:N
    for ii=1:N
        x1 = X1(i,ii);
        x2 = X2(i,ii);
        COSTMIN = 10000;
        UMIN = 10000;
        for jj = U
            Ui = U(jj);
            X_next = A*[x1;x2] + B*Ui;
            C_star = [x1;x2]' * Q * [x1;x2] + Ui' * R * Ui;
            if(C_star < COSTMIN)
                COSTMIN = C_star;
                UMIN = Ui;
            end
        end
        % store UMIN in UOPT(N-k,I)
        u_star(i,ii) = UMIN;
        % store COSMIN in COST(N-k,i)
        J_star(i,ii) = COSTMIN;
    end
end
fprintf('non-vectorized operation complete in %f seconds...\n',toc)

%% ----------------- Vectorized Implementation ---------------------- %%
J_star_V = zeros(N,N);
u_star_V = J_star_V;
JV = zeros(N,N,length(U));
X_next = [(A(1)*X1 + A(3)*X2 + B(1)*Ui);
    (A(2)*X1 + A(4)*X2 + B(2)*Ui)];

C_star = Q(1)*X1.^2 + Q(4)*X2.^2 + R * Ui.^2;


