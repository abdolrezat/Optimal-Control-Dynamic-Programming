% test between different 2d interpolant functions Interp2 and
% griddedinterpolant performances. result: predefined griddedinterpolant is the
% best method.
%% Choose settings
%low dim
N1 = 1000;
N2 = 5000;
N3 = 5000;


[X1,X2] = meshgrid(-100:0.2:100,-100:0.2:100);
J = X1.^2 + 78*X2.^2;
plot3(X1,X2,J)


tic
for ii= 1:N1
    j2 = interp2(X1,X2,J,X1(ii),X2(ii),'linear');
end
toc

X1 = X1';
X2 = X2';

tic
for jj = 1:N2
    F = griddedInterpolant(X1,X2,J,'linear');
    j3 = F(X1(jj),X2(jj));
end
toc

%% new idea
F = griddedInterpolant(X1,X2,J,'linear');
tic
for kk = 1:N3
    J4(kk) = F(X1(kk),X2(kk));
end
toc
hold on

plot3(X1(1:kk), X2(1:kk), J4(1:kk), '*')