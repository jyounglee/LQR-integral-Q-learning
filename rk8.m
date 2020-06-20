function y_n=rk8(fn,y,u,ts,i)
% RK78: Runge-Kutta-Fehlberg method (Runge-Kutta 7th + 8th)
%
%  []=rk8(fn,y,u,ts,i)
%  where, fn is function name recalled
%
% Copyright (C) 2009 Jaeyoung Lee (jyounglee@yonsei.ac.kr)


global rk8_ALPHA;
global rk8_BETA1; global rk8_BETA2; global rk8_BETA3; global rk8_BETA4; global rk8_BETA5;
global rk8_BETA6; global rk8_BETA7; global rk8_BETA8; global rk8_BETA9; global rk8_BETA10;
global rk8_c;

t=(i-1)*ts;
h=ts;

K1=h*feval(fn, y, t, u);
K2=h*feval(fn, y + rk8_BETA1*K1,  t + rk8_ALPHA(1)*h,u);
K3=h*feval(fn, y + rk8_BETA2(1)*K1 + rk8_BETA2(2)*K2,  t+rk8_ALPHA(2)*h,  u);
K4=h*feval(fn, y + rk8_BETA3(1)*K1 + rk8_BETA3(2)*K2 + rk8_BETA3(3)*K3,  t+rk8_ALPHA(3)*h,  u);
K5=h*feval(fn, y + rk8_BETA4(1)*K1 + rk8_BETA4(2)*K2 + rk8_BETA4(3)*K3 + rk8_BETA4(4)*K4,  t+rk8_ALPHA(4)*h,  u);
K6=h*feval(fn, y + rk8_BETA5(1)*K1 + rk8_BETA5(2)*K2 + rk8_BETA5(3)*K3 + rk8_BETA5(4)*K4 + rk8_BETA5(5)*K5,  t+rk8_ALPHA(5)*h,  u);
K7=h*feval(fn, y + rk8_BETA6(1)*K1 + rk8_BETA6(2)*K2 + rk8_BETA6(3)*K3 + rk8_BETA6(4)*K4 + rk8_BETA6(5)*K5 + rk8_BETA6(6)*K6,  t+rk8_ALPHA(6)*h,  u);
K8=h*feval(fn, y + rk8_BETA7(1)*K1 + rk8_BETA7(2)*K2 + rk8_BETA7(3)*K3 + rk8_BETA7(4)*K4 + rk8_BETA7(5)*K5 + rk8_BETA7(6)*K6 + rk8_BETA7(7)*K7,  t+rk8_ALPHA(7)*h,  u);
K9=h*feval(fn, y + rk8_BETA8(1)*K1 + rk8_BETA8(2)*K2 + rk8_BETA8(3)*K3 + rk8_BETA8(4)*K4 + rk8_BETA8(5)*K5 + rk8_BETA8(6)*K6 + rk8_BETA8(7)*K7 + rk8_BETA8(8)*K8,  t+rk8_ALPHA(8)*h,  u);
K10=h*feval(fn, y + rk8_BETA9(1)*K1 + rk8_BETA9(2)*K2 + rk8_BETA9(3)*K3 + rk8_BETA9(4)*K4 + rk8_BETA9(5)*K5 + rk8_BETA9(6)*K6 + rk8_BETA9(7)*K7 + rk8_BETA9(8)*K8 + rk8_BETA9(9)*K9,  t+rk8_ALPHA(9)*h,  u);
K11=h*feval(fn, y + rk8_BETA10(1)*K1 + rk8_BETA10(2)*K2 + rk8_BETA10(3)*K3 + rk8_BETA10(4)*K4 + rk8_BETA10(5)*K5 + rk8_BETA10(6)*K6 + rk8_BETA10(7)*K7 + rk8_BETA10(8)*K8 + rk8_BETA10(9)*K9 + rk8_BETA10(10)*K10,  t+rk8_ALPHA(10)*h,  u);

y_n = y + rk8_c*[K1;K8;K9;K10;K11];
