function rk8_init
% rk8_init: Initialization function for rk8
%
% Call this function before you use rk8(.)
%
% Copyright (C) 2009 Jaeyoung Lee (jyounglee@yonsei.ac.kr)

global rk8_ALPHA;
global rk8_BETA1; global rk8_BETA2; global rk8_BETA3; global rk8_BETA4; global rk8_BETA5;
global rk8_BETA6; global rk8_BETA7; global rk8_BETA8; global rk8_BETA9; global rk8_BETA10;
global rk8_c;

s21 = sqrt(21);
rk8_ALPHA = [1/2 1/2 (7+s21)/14 (7+s21)/14 1/2 (7-s21)/14 (7-s21)/14 1/2 (7+s21)/14 1];
rk8_BETA1 = 1/2;
rk8_BETA2 = [1/4 1/4];                               
rk8_BETA3 = [1/7 (-7-3*s21)/98 (21+5*s21)/49]; 
rk8_BETA4 = [(11+s21)/84 0 (18+4*s21)/63 (21-s21)/252]; 
rk8_BETA5 = [(5+s21)/48 0 (9+s21)/36  (-231+14*s21)/360 (63-7*s21)/80]; 
rk8_BETA6 = [(10-s21)/42 0 (-432+92*s21)/315 (633-145*s21)/90   (-504+115*s21)/70 (63-13*s21)/35 ];
rk8_BETA7 = [1/14 0 0 0 (14-3*s21)/126 (13-3*s21)/63 1/9];
rk8_BETA8 = [1/32 0 0 0 (91-21*s21)/576  11/72  (-385-75*s21)/1152  (63+13*s21)/128];
rk8_BETA9 = [1/14  0 0 0 1/9 (-733-147*s21)/2205 (515+111*s21)/504  (-51-11*s21)/56  (132+28*s21)/245];
rk8_BETA10 = [0 0 0 0 (-42+7*s21)/18  (-18+28*s21)/45  (-273-53*s21)/72 (301+53*s21)/72  (28-28*s21)/45  (49-7*s21)/18];
rk8_c = [9 49 64 49 9]./180;  % y[i+1] = y[i] + ( 9*k1 + 49*k8 + 64*k9 + 49*k10 + 9*k11 ) / 180
