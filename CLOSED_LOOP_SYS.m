function x_dot = CLOSED_LOOP_SYS(x_in, u, t)
% Calculate the derivative of the system and learning states
%
% Copyright (C) 2009 Jaeyoung Lee (jyounglee@yonsei.ac.kr)

global A;
global B;
global Q;
global R;
global ws;

w = ws;

x_dot = zeros(9,1);
x_dot(1:4) = A*x_in(1:4)' + B*(u+w);
x_dot(5) =  x_in(1:4)*Q*x_in(1:4)' + u*R*u';
x_dot(6) =  x_in(1)*w;
x_dot(7) =  x_in(2)*w;
x_dot(8) =  x_in(3)*w;
x_dot(9) =  x_in(4)*w;

x_dot = x_dot';
