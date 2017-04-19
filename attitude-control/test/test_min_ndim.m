u1 = [3 -2 1];
u2 = [1 2 -2];

x1 = [1 2 3]
[X1,U1,U2] = ndgrid(x1,u1,u2);
J = X1.*U1.*U2;

[val, id2] = min(J,[],3)
[val, id1] = min(val,[],2)

u2(id2(id1))
u1(id1)
x.*u2(id2(id1)).*u1(id1)

