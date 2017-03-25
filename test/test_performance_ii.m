% tests to assess how much performance is improved if a single line of
% allocation e.g. x1 = X1(i,ii) executed in the outer loop. of course the
% results will differ cause the logics are not the same

[X1, X2] = ndgrid(1:100,1:100);
U = 1:1000;

J = zeros(100,100,1000);

for k = 1:10
    tic
    for i=1:100
        for ii=1:100
            x1 = X1(i,ii);
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