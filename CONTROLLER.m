function u = CONTROLLER(x_in)
% Calculate the feedback control input.
%
% Copyright (C) 2009 Jaeyoung Lee (jyounglee@yonsei.ac.kr)

 global K_x; 
 u = K_x(1)*x_in(1) + K_x(2)*x_in(2) + K_x(3)*x_in(3) + K_x(4)*x_in(4);
