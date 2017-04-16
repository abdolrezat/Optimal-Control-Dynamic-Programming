%% test gridded interpolant see if it works as desired (works)
obj = Dynamic_Solver;
s_r = linspace(obj.x_min,obj.x_max,obj.dx);
[obj.X1_mesh, obj.X2_mesh] = ndgrid(s_r, s_r);
U_mesh = linspace(obj.u_max, obj.u_min,obj.du);
obj.J_star = zeros([size(obj.X1_mesh),obj.N]);
obj.u_star = obj.J_star;
k  = 1;

%make a surface
obj.J_star(:,:,obj.N-k+1) = 2.* obj.X1_mesh .* obj.X2_mesh + obj.X2_mesh
obj.J_star(:,:,end) %display surface Z
F = griddedInterpolant(obj.X1_mesh, obj.X2_mesh,...
obj.J_star(:,:,obj.N-k+1),'linear');

plot3(obj.X1_mesh, obj.X2_mesh, obj.J_star(:,:,end))
hold on
plot3(0,0,F(0,0),'*') %plot a point
p = @(a,b)(plot3(a,b,F(a,b),'*r')); %function to plot points
p([7 5 3 3 1 -4 -5],[10 8 5 3 -1 -2 -5]) %plot random points see if they fit the surface
%it works

%test with only grid vectors {s_r,s_r}
F22 = griddedInterpolant({s_r,s_r}, obj.J_star(:,:,end));

%test memory
J = obj.J_star(:,:,end);
whos F22 J

%test one Y(1xN) vector 
%FY = griddedInterpolant({s_r,s_r}, J(:)')
%FXY = griddedInterpolant(obj.X1_mesh(:), obj.X2_mesh(:), J(:),'linear')

%test performance of 'nearest' vs. 'linear'

FN = griddedInterpolant(obj.X1_mesh, obj.X2_mesh,...
obj.J_star(:,:,obj.N-k+1),'nearest');

pN = @(a,b)(plot3(a,b,FN(a,b),'*b')); %function to plot points
pN([7 5 3 3 1 -4 -5],[10 8 5 3 -1 -2 -5]) %plot random points see if they fit the surface


tic
X = J;
for ii=1:10000
    X = F(obj.X1_mesh, obj.X2_mesh);
end
toc
tic
X = J;
for ii=1:10000
    X = FN(obj.X1_mesh, obj.X2_mesh);
end
toc

