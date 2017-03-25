%tests the performance difference if X1(i,ii) and other matrices are
%implemented directly in loops instead of allocation to another variable
%e.g. x1 = X1(i,ii);
% Result: time is about 50% reduced when the matrix element is 
[X1, X2] = ndgrid(1:100,1:100);
U = 1:1000;

J = zeros(100,100,1000);

for k = 1:10
    tic
    for i=1:100
        for ii=1:100
            COSTMIN = Inf;
            for j = 1:1000
                C_star = X1(i,ii)*2 + 2*X2(i,ii)^2 + 300*U(j);
                if C_star < COSTMIN
                    COSTMIN = C_star;
                end
            end
            J(i,ii,k) = COSTMIN;
        end
    end
    fprintf('test 1: step %d - %f seconds\n', k, toc)
end

%% ------------------------------------ Second Test ----------------------------- %%


[X1, X2] = ndgrid(1:100,1:100);
U = 1:1000;

J = zeros(100,100,1000);

for k = 1:10
    tic
    for i=1:100
        x1 = X1(i,ii);
        for ii=1:100
            x2 = X2(i,ii);
            COSTMIN = Inf;
            for j = 1:1000
                u = U(j);
                C_star = x1*2 + 2*x2^2 + 300*u;
                if C_star < COSTMIN
                    COSTMIN = C_star;
                end
            end
            J(i,ii,k) = COSTMIN;
        end
    end
    fprintf('test 2: step %d - %f seconds\n', k, toc)
end