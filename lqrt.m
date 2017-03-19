function [K,S]=lqrt(F,Q,R,A,B,t,tf)

%  	[K,S]=lqrt(F,Q,R,A,B,t,tf)
%	Linear quadratic regulator design for continuous systems.
%	Calculates the optimal feedback gain (time varing)
%	matrix K(t) such that the feedback law  u = -K(t)x  minimizes the cost
%	function:
%
%		J =x(tf)'Fx(tf) + Integral {x'Qx + u'Ru} dt
%
%	subject to the constraint equation: 
%		.
%		x = Ax + Bu 
%                
%	Also returned is S(t), the solution to the associated Riccati equation 
%	.		     -1
%	S + SA + A'S - SB * R * B'S + Q = 0    
%
%	See also: LQR, LQRY, LQR2, and REG.

%       Written by Giampiero Campa , 6-11-94.

[n,n]=size(A);

H=[A -B*inv(R)*B'; -Q -A'];

fi=expm(H*(tf-t));
fi11=fi(1:n,1:n);
fi12=fi(1:n,n+1:2*n);
fi21=fi(n+1:2*n,1:n);
fi22=fi(n+1:2*n,n+1:2*n);

S=inv(fi22-F*fi12)*(fi21-F*fi11);
K=inv(R)*B'*S;
