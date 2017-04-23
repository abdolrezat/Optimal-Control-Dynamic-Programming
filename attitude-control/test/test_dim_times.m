%%
x = [1,2]
u1 = [3 0 -2]'
u2 = reshape([-4 2 4],[1 1 3])

J1 = x + u1 + 3*u2

%%
[X,U1,U2] = ndgrid([3 0 -2],[1,2],[-4 2 4])
J2 = X + U1 + 3*U2
%%

J3_p = u1 + 3*u2

%%

yaw1 = [0,1,2]
pitch1 = [1,1.5,2]
roll1 = [-1,-2,5]

%% n-D gridded implementation / uses memory to store grid
[X1,X2,X5,X6,X7] = ndgrid(x,x,yaw1,pitch1,roll1);
size_Q = size(X5);
angles = [X5(:) X6(:) X7(:)];

cang = cos( angles/2 );
sang = sin( angles/2 );
%start test
tic
X5 = reshape( cang(:,1).*sang(:,2).*cang(:,3) + sang(:,1).*cang(:,2).*sang(:,3), size_Q);
X6 = reshape( cang(:,1).*cang(:,2).*sang(:,3) - sang(:,1).*sang(:,2).*cang(:,3), size_Q);
X7 = reshape( cang(:,1).*cang(:,2).*cang(:,3) + sang(:,1).*sang(:,2).*sang(:,3), size_Q);

JQ1 = X1.^2 + 2*X2.^2 + X5*2 + X6*3 + X7*5;
toc
%end test

%% vector implementation / saves memory
cang_x5 = reshape( cos( yaw1/2 ), [1 1 3 1 1]);
sang_x5 = reshape( sin( yaw1/2 ), [1 1 3 1 1]);

cang_x6 = reshape( cos( pitch1/2 ), [1 1 1 3 1]);
sang_x6 = reshape( sin( pitch1/2 ), [1 1 1 3 1]);

cang_x7 = reshape( cos( roll1/2 ), [1 1 1 1 3]);
sang_x7 = reshape( sin( roll1/2 ), [1 1 1 1 3]);
X1V = reshape(x,[2 1 1 1 1]);
X2V = reshape(x,[1 2 1 1 1]);
%start test
tic
X5V = cang_x5.*sang_x6.*cang_x7 + sang_x5.*cang_x6.*sang_x7;
X6V = cang_x5.*cang_x6.*sang_x7 - sang_x5.*sang_x6.*cang_x7;
X7V = cang_x5.*cang_x6.*cang_x7 + sang_x5.*sang_x6.*sang_x7;
JQ2 = X1V.^2 + 2*X2V.^2 + X5V*2 + X6V*3 + X7V*5;
toc
%end test
%% Direct vector implementation / much less memory and time required
%start test
tic
JQ3 = X1V.^2 + 2*X2V.^2 + ...
    (cang_x5.*sang_x6.*cang_x7 + sang_x5.*cang_x6.*sang_x7)*2 + ...
    (cang_x5.*cang_x6.*sang_x7 - sang_x5.*sang_x6.*cang_x7)*3 + ...
    (cang_x5.*cang_x6.*cang_x7 + sang_x5.*sang_x6.*sang_x7)*5;
toc
%end test
