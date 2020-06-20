
% Main program to reproduce integral Q-learning, Fig 1(a),(d),
% in the publication:
% 
% Lee, J.Y., Park, J.B., and Choi, Y.H.,
% Integral Q-learning and explorized policy iteration
% for adaptive optimal control of continuous-time linear systems,
% Automatica 11(48), 2850--2859, 2012.
%
% Copyright (C) 2009 Jaeyoung Lee (jyounglee@yonsei.ac.kr)

clc;
clear all;
close all;

format

t0 = 0; tf = 7.8;

ts = 0.001; % Inner Sampling Time
Ts = 0.01;  % Outer Sampling Time
Tadp = 0.05; % Iteration Period

InnerIter = fix(Ts/ts);
OuterIter = fix(Tadp/Ts);
ADPIter = fix((tf-t0)/Tadp);

max_update = 9;

t = 0:Ts:tf-Ts;
t_temp = 0:ts:Ts;

rk8_init;

global R;
global Q;
global A;
global B;
global K_x;
global HNOW;
global ws;

A = [-0.0665  11.5        0         0 ; ...
    0     -2.5       2.5        0 ; ...
    -9.5       0   -13.736  -13.736 ; ...
    0.6       0        0         0];
B = [0 0 13.736 0]';
Q = eye(4);
R = 1;

N = 14;
N_min = 14;

[K, P] = lqr(A,B,Q,R,0);
h = zeros(14,1);
HNOW = h;
K_x = [0 0 0 0];

V = 0;
W = [0 0 0 0];

x0 = [0 0 0 0]';
x_temp = x0';
u_temp = K_x*x0;

x = [x0 zeros(4,OuterIter*ADPIter-1)];
u = [u_temp zeros(1,OuterIter*ADPIter-1)];

X = [x0(1)^2 x0(1)*x0(2) x0(1)*x0(3) x0(1)*x0(4) ...
    x0(2)^2  x0(2)*x0(3) x0(2)*x0(4) ...
    x0(3)^2  x0(3)*x0(4) ...
    x0(4)^2]';
dX = zeros(10,1);
dZ = zeros(14,1);

Y = [];

l = 1;  pSize = 1;  r = 0;  condn = [];
learning_points = [t0 ; x0 ; u_temp];

x_m = 10^(-3); w_m = 10^(-3);

for i = 1:ADPIter
    ws = w_m*(2 * randn - 1);
    for j = 1:OuterIter
        for k = 1:InnerIter
            t_now = (((i-1)*OuterIter + j-1)*InnerIter + k-1)*ts;
            [x_temp, V, W] = CalculateNext(x_temp, V, W, t_now, ts);
            u_temp = CONTROLLER(x_temp);
        end
        OUTERINDEX = (i-1)*OuterIter+j;
        x(:,OUTERINDEX) = x_temp;
        u(OUTERINDEX) = CONTROLLER(x_temp);
    end
    
    % Data Aquisition for Performing Least Squares %%%%%%%%%%%%%%%%%%%%%%%
    dX(:,l) = X - [x_temp(1)^2 x_temp(1)*x_temp(2)   x_temp(1)*x_temp(3)   x_temp(1)*x_temp(4) ...
        x_temp(2)^2 x_temp(2)*x_temp(3)   x_temp(2)*x_temp(4) ...
        x_temp(3)^2 x_temp(3)*x_temp(4) ...
        x_temp(4)^2]';
    X = X - dX(:,l);
    dZ(:,l) = [dX(:,l)' 2*W(1) 2*W(2) 2*W(3) 2*W(4)]';
    Y(l) = V;
    V = 0;
    W = [0 0 0 0];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    r = rank(dZ(:,1:l));
    
    if(r == N_min && l >= N && pSize <= max_update)
        [Kron_U Kron_S Kron_V] = svd(dZ(:,1:l)');
        condn = [condn, Kron_S(1,1)/Kron_S(r,r)];
        disp(['Condition Number of the Matrix: ', num2str(Kron_S(1,1)/Kron_S(r,r)), ' (iter: ', num2str(i), ')'])
        [Ssrow Sscol] = size(Kron_S);
        Inv_Kron_S = inv(Kron_S(1:r,1:r));
        Inv_Kron_S = [Inv_Kron_S zeros(r,  Ssrow-r)];
        Inv_Kron_S = [Inv_Kron_S ; zeros(Sscol-r, Ssrow)];
        h(:,pSize+1) = Kron_V*Inv_Kron_S*Kron_U'*Y';
        pSize = pSize+1;
        
        % Information for ploting the graphs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        learning_points(:,pSize) = [t(OUTERINDEX); x(:,OUTERINDEX); u(OUTERINDEX)];
        
        % Controller & Value  Update %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        HNOW = h(:,pSize);
        PU = [HNOW(1) HNOW(2) HNOW(3)  HNOW(4); ...
            0   HNOW(5) HNOW(6)  HNOW(7); ...
            0        0  HNOW(8)  HNOW(9); ...
            0        0        0  HNOW(10)];
        PU = (PU + PU')/2;
        
        H21 = [HNOW(11) HNOW(12) HNOW(13) HNOW(14)];
        K_x = -inv(R)*H21;
        
        % Initialization for the next update %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Y = [];
        r = 0;        l = 0;
        
        if pSize > max_update
            w_m = 0;
        else
            M = Q + K_x'* R * K_x;
            Pi = [M H21' ; H21 2*R];
            C = norm(H21) / min(eig(M));
            D = norm([H21  R]) / min(eig(Pi));
            w_m = x_m / min(C, D) / 2;
        end
    end
    l = l+1;
end

LineW = 1.5;

%% Ploting the trj's of states %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(t,x(1,:), 'LineWidth', 2, 'Color', [0 0.8 0.8]);
hold on
plot(t,x(2,:), '--', 'LineWidth', 2, 'Color', [0 0.7 0]);
hold on
plot(t,x(3,:), ':', 'LineWidth', 2, 'Color', [1 0 0]);
hold on
plot(t,x(4,:), '-.', 'LineWidth', 2, 'Color', [0 0 0.7]);

n = pSize;
while (n>1)
    plot(learning_points(1, n), learning_points(2, n)','.', 'MarkerSize', 20, 'Color', [0 0.8 0.8]);
    hold on
    plot(learning_points(1, n), learning_points(3, n)','.', 'MarkerSize', 20, 'Color', [0 0.7 0]);
    hold on
    plot(learning_points(1, n), learning_points(4, n)','.', 'MarkerSize', 20, 'Color', [1 0 0]);
    hold on
    plot(learning_points(1, n), learning_points(5, n)','.', 'MarkerSize', 20, 'Color', [0 0 0.7]);
    hold on
    n = n-1;
end

% for n=1:ADPIter
%      tPoint = n*OuterIter - 1;
%      plot(t(tPoint), x(1,tPoint),'.', 'MarkerSize', 10, 'Color','red');
%      hold on
%      plot(t(tPoint), x(2,tPoint),'.', 'MarkerSize', 10, 'Color', [0 0 1]);
%      hold on
%      plot(t(tPoint), x(3,tPoint),'.', 'MarkerSize', 10, 'Color', [0 0.5 0]);
%      hold on
%      plot(t(tPoint), x(4,tPoint),'.', 'MarkerSize', 10, 'Color', [0 0.8 0.8]);
%      hold on
%      n = n-1;
% end

title('State Trajectories');
xlabel('Time(sec)');
legend('x_1', 'x_2', 'x_3', 'x_4');
xlim([0 tf]);

%% Ploting critic parameters variations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
l = length(learning_points(1, :));
plot(learning_points(1, 1:pSize), h(1,1:pSize)' ,'bd--', 'LineWidth',LineW, 'MarkerSize', 9);
hold on
plot(learning_points(1, 1:pSize), h(8,1:pSize)' ,'.--', 'LineWidth',LineW,  'MarkerSize', 20, 'Color', [0 0.8 0.8]);
hold on
plot(learning_points(1, 1:pSize), h(6,1:pSize)'./2 ,'gx--', 'LineWidth',LineW,  'MarkerSize', 9);
hold on
plot(learning_points(1, 1:pSize), h(7,1:pSize)'./2 ,'rs--', 'LineWidth',LineW,  'MarkerSize', 9);


legend('\it{P}_{11}', '\it{P}_{33}', '\it{P}_{23}', '\it{P}_{42}');
xlabel('Time(sec)');
title('Evolution of the Parameters of \it P');
xlim([0 6.5]);
ylim([-0.01 3.6]);

%% Ploting actor parameters variations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
l = length(learning_points(1, :));
plot(learning_points(1, 1:pSize), h(11,1:pSize)' ,'gx--', 'LineWidth',LineW,  'MarkerSize', 9);
hold on
plot(learning_points(1, 1:pSize), h(12,1:pSize)' ,'rs--', 'LineWidth',LineW,  'MarkerSize', 9);
hold on
plot(learning_points(1, 1:pSize), h(13,1:pSize)' ,'bd--', 'LineWidth',LineW, 'MarkerSize', 9);
hold on
plot(learning_points(1, 1:pSize), h(14,1:pSize)' ,'.--', 'LineWidth',LineW,  'MarkerSize', 20, 'Color', [0 0.8 0.8]);
hold on

legend('\it{K}_{1}', '\it{K}_{2}', '\it{K}_{3}', '\it{K}_{4}');
xlabel('Time(sec)');
title('Evolution of the Parameters of \it K');
xlim([0 6.5]);
ylim([-0.1 41]);

%% Ploting control input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
n = pSize;
plot(t((learning_points(1, n)/Ts+2):length(t)), u((learning_points(1, n)/Ts+2):length(u)), 'k', 'LineWidth', LineW);
hold on
plot(t(int32(learning_points(1, n)/Ts)+1), u(int32(learning_points(1, n)/Ts)+1),'k.:', 'MarkerSize', 15);
hold on

while (n>2)
    plot(t((learning_points(1, n-1)/Ts+2):learning_points(1, n)/Ts), u((learning_points(1, n-1)/Ts+2):learning_points(1, n)/Ts), 'k', 'LineWidth',LineW);
    hold on
    plot(t(int32(learning_points(1, n-1)/Ts)+2), u(int32(learning_points(1, n-1)/Ts)+2),'k.:', 'MarkerSize', 15);
    hold on
    n = n -1;
end
hold on
plot(t((learning_points(1, n-1)/Ts+2):learning_points(1, n)/Ts), u((learning_points(1, n-1)/Ts+2):learning_points(1, n)/Ts), 'k', 'LineWidth',LineW);


xlabel('Time(sec)');
ylabel('Control input u');
title('The trajectory of the input variable u');

