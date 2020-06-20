function [x_n, V_n, W_n] = CalculateNext(x,V,W,t,ts)
% Calculate the next system and learning states
% using Runge-Kutta-Fehlberg method (Runge-Kutta 7th + 8th)
%
% Copyright (C) 2009 Jaeyoung Lee (jyounglee@yonsei.ac.kr)

global rk8_ALPHA;
global rk8_BETA1; global rk8_BETA2; global rk8_BETA3; global rk8_BETA4; global rk8_BETA5;
global rk8_BETA6; global rk8_BETA7; global rk8_BETA8; global rk8_BETA9; global rk8_BETA10;
global rk8_c;

h=ts;
X = [x V W];

K1=h*feval('CLOSED_LOOP_SYS', X, feval('CONTROLLER', X),t);
X1 = X + rk8_BETA1*K1;

K2=h*feval('CLOSED_LOOP_SYS', X1, feval('CONTROLLER', X1),t + rk8_ALPHA(1)*h);
X2 = X + rk8_BETA2(1)*K1 + rk8_BETA2(2)*K2;

K3=h*feval('CLOSED_LOOP_SYS', X2, feval('CONTROLLER', X2), t + rk8_ALPHA(2)*h);
X3 = X + rk8_BETA3(1)*K1 + rk8_BETA3(2)*K2 + rk8_BETA3(3)*K3;

K4=h*feval('CLOSED_LOOP_SYS', X3, feval('CONTROLLER', X3), t + rk8_ALPHA(3)*h);
X4 = X + rk8_BETA4(1)*K1 + rk8_BETA4(2)*K2 + rk8_BETA4(3)*K3 + rk8_BETA4(4)*K4;

K5=h*feval('CLOSED_LOOP_SYS', X4, feval('CONTROLLER', X4), t + rk8_ALPHA(4)*h);
X5 = X + rk8_BETA5(1)*K1 + rk8_BETA5(2)*K2 + rk8_BETA5(3)*K3 + rk8_BETA5(4)*K4 + rk8_BETA5(5)*K5;

K6=h*feval('CLOSED_LOOP_SYS', X5, feval('CONTROLLER', X5), t + rk8_ALPHA(5)*h);
X6 = X + rk8_BETA6(1)*K1 + rk8_BETA6(2)*K2 + rk8_BETA6(3)*K3 + rk8_BETA6(4)*K4 + rk8_BETA6(5)*K5 + rk8_BETA6(6)*K6;

K7=h*feval('CLOSED_LOOP_SYS', X6, feval('CONTROLLER', X6), t + rk8_ALPHA(6)*h);
X7 = X + rk8_BETA7(1)*K1 + rk8_BETA7(2)*K2 + rk8_BETA7(3)*K3 + rk8_BETA7(4)*K4 + rk8_BETA7(5)*K5 + rk8_BETA7(6)*K6 + rk8_BETA7(7)*K7;

K8=h*feval('CLOSED_LOOP_SYS', X7, feval('CONTROLLER', X7), t + rk8_ALPHA(7)*h);
X8 = X + rk8_BETA8(1)*K1 + rk8_BETA8(2)*K2 + rk8_BETA8(3)*K3 + rk8_BETA8(4)*K4 + rk8_BETA8(5)*K5 + rk8_BETA8(6)*K6 + rk8_BETA8(7)*K7 + rk8_BETA8(8)*K8;

K9=h*feval('CLOSED_LOOP_SYS', X8, feval('CONTROLLER', X8), t + rk8_ALPHA(8)*h);
X9 = X + rk8_BETA9(1)*K1 + rk8_BETA9(2)*K2 + rk8_BETA9(3)*K3 + rk8_BETA9(4)*K4 + rk8_BETA9(5)*K5 + rk8_BETA9(6)*K6 + rk8_BETA9(7)*K7 + rk8_BETA9(8)*K8 + rk8_BETA9(9)*K9;

K10=h*feval('CLOSED_LOOP_SYS', X9, feval('CONTROLLER', X9), t + rk8_ALPHA(9)*h);
X10 = X + rk8_BETA10(1)*K1 + rk8_BETA10(2)*K2 + rk8_BETA10(3)*K3 + rk8_BETA10(4)*K4 + rk8_BETA10(5)*K5 + rk8_BETA10(6)*K6 + rk8_BETA10(7)*K7 + rk8_BETA10(8)*K8 + rk8_BETA10(9)*K9 + rk8_BETA10(10)*K10;

K11=h*feval('CLOSED_LOOP_SYS', X10, feval('CONTROLLER', X10), t + rk8_ALPHA(10)*h);

X_n = X + rk8_c*[K1;K8;K9;K10;K11];
x_n = X_n(1:4);
V_n = X_n(5);
W_n = X_n(6:9);
